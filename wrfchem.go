/*
Copyright © 2013 the InMAP authors.
This file is part of InMAP.

InMAP is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

InMAP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with InMAP.  If not, see <http://www.gnu.org/licenses/>.
*/

package inmap

import (
	"fmt"
	"math"
	"time"

	"github.com/ctessum/atmos/seinfeld"
	"github.com/ctessum/atmos/wesely1989"

	"github.com/ctessum/sparse"
)

// WRF variables currently used:
/* hc5,hc8,olt,oli,tol,xyl,csl,cvasoa1,cvasoa2,cvasoa3,cvasoa4,iso,api,sesq,lim,
cvbsoa1,cvbsoa2,cvbsoa3,cvbsoa4,asoa1i,asoa1j,asoa2i,asoa2j,asoa3i,asoa3j,asoa4i,
asoa4j,bsoa1i,bsoa1j,bsoa2i,bsoa2j,bsoa3i,bsoa3j,bsoa4i,bsoa4j,no,no2,no3ai,no3aj,
so2,sulf,so4ai,so4aj,nh3,nh4ai,nh4aj,PM2_5_DRY,U,V,W,PBLH,PH,PHB,HFX,UST,PBLH,T,
PB,P,ho,h2o2,LU_INDEX,QRAIN,CLDFRA,QCLOUD,ALT,SWDOWN,GLW */

const wrfFormat = "2006-01-02_15_04_05"

// WRFChem is an InMAP preprocessor for WRF-Chem output.
type WRFChem struct {
	aVOC, bVOC, aSOA, bSOA, nox, no, no2, pNO, sox, pS, nh3, pNH, totalPM25 map[string]float64

	start, end time.Time

	wrfOut string

	recordDelta, fileDelta time.Duration

	msgChan chan string
}

// NewWRFChem initializes a WRF-Chem preprocessor from the given
// configuration information.
// WRFOut is the location of WRF-Chem output files.
// [DATE] should be used as a wild card for the simulation date.
// startDate and endDate are the dates of the beginning and end of the
// simulation, respectively, in the format "YYYYMMDD".
// If msgChan is not nil, status messages will be sent to it.
func NewWRFChem(WRFOut, startDate, endDate string, msgChan chan string) (*WRFChem, error) {
	w := WRFChem{
		// These maps contain the WRF-Chem variables that make
		// up the chemical species groups, as well as the
		// multiplication factors required to convert concentrations
		// to mass fractions [μg/kg dry air].

		// RACM VOC species and molecular weights (g/mol);
		// Only includes anthropogenic precursors to SOA from
		// anthropogenic (aSOA) and biogenic (bSOA) sources as
		// in Ahmadov et al. (2012)
		// Assume condensable vapor from SOA has molar mass of 70
		aVOC: map[string]float64{"SOAALK": ppmvToUgKg(112), "TOLUENE": ppmvToUgKg(92.14),
			"MXYL": ppmvToUgKg(106.17), "OXYL": ppmvToUgKg(106.17), "PXYL": ppmvToUgKg(106.17),
			"ARO1": ppmvToUgKg(95.16), "ARO2MN": ppmvToUgKg(118.72),
			"NAPHTHAL": ppmvToUgKg(128.20), "TMBENZ124": ppmvToUgKg(120.19),
			"ISOPRENE": ppmvToUgKg(68.12), "APIN": ppmvToUgKg(136.23), "SESQ": ppmvToUgKg(204.35),
			"TERP": ppmvToUgKg(136.24)},
		bVOC: map[string]float64{
			"ISOPRENE": ppmvToUgKg(68.12), "APIN": ppmvToUgKg(136.23), "SESQ": ppmvToUgKg(204.35),
			"TERP": ppmvToUgKg(136.24)},
		// VBS SOA species (anthropogenic only) [μg/m3].
		aSOA: map[string]float64{"AAVB1J": 1, "AAVB2J": 1, "AAVB3J": 1,
			"AAVB4J": 1, "AOLGAJ": 1, "APCSOJ": 1, "AISO1J": 1, "AISO2J": 1, "AISO3J": 1,
			"ASQTJ": 1, "AOLGBJ": 1, "AMT1J": 1, "AMT2J": 1, "AMT3J": 1,
			"AMT4J": 1, "AMT5J": 1, "AMT6J": 1, "AIETETJ": 1, "AIEOSJ": 1,
			"ADIMJ": 1, "AIMGAJ": 1, "AIMOSJ": 1, "AMTNO3J": 1, "AISOPNNJ": 1,
			"AMTHYDJ": 1, "AGLYJ": 1},
		// VBS SOA species (biogenic only) [μg/m3].
		bSOA: map[string]float64{"AISO1J": 1, "AISO2J": 1, "AISO3J": 1,
			"ASQTJ": 1, "AOLGBJ": 1, "AMT1J": 1, "AMT2J": 1, "AMT3J": 1,
			"AMT4J": 1, "AMT5J": 1, "AMT6J": 1, "AIETETJ": 1, "AIEOSJ": 1,
			"ADIMJ": 1, "AIMGAJ": 1, "AIMOSJ": 1, "AMTNO3J": 1, "AISOPNNJ": 1,
			"AMTHYDJ": 1, "AGLYJ": 1},
		// NOx is RACM NOx species. We are only interested in the mass
		// of Nitrogen, rather than the mass of the whole molecule, so
		// we use the molecular weight of Nitrogen.
		nox: map[string]float64{"NO": ppmvToUgKg(mwN), "NO2": ppmvToUgKg(mwN)},
		// pNO is the Nitrogen fraction of MADE particulate
		// NO species [μg/m3].
		pNO: map[string]float64{"ANO3I": mwN / mwNO3, "ANO3J": mwN / mwNO3},
		// SOx is the RACM SOx species. We are only interested in the mass
		// of Sulfur, rather than the mass of the whole molecule, so
		// we use the molecular weight of Sulfur.
		sox: map[string]float64{"SO2": ppmvToUgKg(mwS), "SULF": ppmvToUgKg(mwS)},
		// pS is the Sulfur fraction of the MADE particulate
		// Sulfur species [μg/m3].
		pS: map[string]float64{"ASO4I": mwS / mwSO4, "ASO4J": mwS / mwSO4},
		// NH3 is ammonia. We are only interested in the mass
		// of Nitrogen, rather than the mass of the whole molecule, so
		// we use the molecular weight of Nitrogen.
		nh3: map[string]float64{"NH3": ppmvToUgKg(mwN)},
		// pNH is the Nitrogen fraction of the MADE particulate
		// ammonia species [μg/m3].
		pNH: map[string]float64{"ANH4I": mwN / mwNH4, "ANH4J": mwN / mwNH4},
		// totalPM25 is total mass of PM2.5  [μg/m3].
		totalPM25: map[string]float64{"ATOTIJ": 1},

		//totalPM25: map[string]float64{"AAVB1J": 1, "AAVB2J": 1, "AAVB3J": 1,
		//	"AAVB4J": 1, "AOLGAJ": 1, "AISO1J": 1, "AISO2J": 1, "AISO3J": 1,
		//	"ASQTJ": 1, "AOLGBJ": 1, "AMT1J": 1, "AMT2J": 1, "AMT3J": 1,
		//	"AMT4J": 1, "AMT5J": 1, "AMT6J": 1, "AIETETJ": 1, "AIEOSJ": 1,
		//	"ADIMJ": 1, "AIMGAJ": 1, "AIMOSJ": 1, "AMTNO3J": 1, "AISOPNNJ": 1,
		//	"AMTHYDJ": 1, "AGLYJ": 1, "ANO3I": 1, "ANO3J": 1,
		//	"ASO4I": 1, "ASO4J": 1, "ANH4I": 1, "ANH4J": 1,
		//	"ANAI": 1, "ANAJ": 1, "ACLI": 1, "ACLJ": 1, "AECI": 1, "AECJ": 1, "AOTHRI": 1, "AOTHRJ": 1, "AORGCJ": 1,
		//	"AFEJ": 1, "ASIJ": 1, "ATIJ": 1, "ACAJ": 1, "AMGJ": 1, "AMNJ": 1, "AALJ": 1, "AKJ": 1,
		//	"ALVPO1I": 1, "ALVPO1J": 1, "ASVPO1I": 1, "ASVPO1J": 1, "ASVPO2I": 1, "ASVPO2J": 1, "ASVPO3J": 1,
		//	"AIVPO1J": 1, "ALVOO1I": 1, "ALVOO1J": 1, "ALVOO2I": 1, "ALVOO2J": 1, "ASVOO1I": 1, "ASVOO1J": 1,
		//	"ASVOO2I": 1, "ASVOO2J": 1, "ASVOO3J": 1}, //APCSOJ-SOA, APOC-POA, APNCOM-POA
		//ACORS, ASOIL, ASEACAT

		wrfOut:  WRFOut,
		msgChan: msgChan,
	}

	var err error
	w.start, err = time.Parse(inDateFormat, startDate)
	if err != nil {
		return nil, fmt.Errorf("inmap: WRF-Chem preprocessor start time: %v", err)
	}
	w.end, err = time.Parse(inDateFormat, endDate)
	if err != nil {
		return nil, fmt.Errorf("inmap: WRF-Chem preprocessor end time: %v", err)
	}

	if !w.end.After(w.start) {
		if err != nil {
			return nil, fmt.Errorf("inmap: WRF-Chem preprocessor end time %v is not after start time %v", w.end, w.start)
		}
	}

	w.recordDelta, err = time.ParseDuration("1h")
	if err != nil {
		return nil, fmt.Errorf("inmap: WRF-Chem preprocessor recordDelta: %v", err)
	}
	w.fileDelta, err = time.ParseDuration("24h")
	if err != nil {
		return nil, fmt.Errorf("inmap: WRF-Chem preprocessor fileDelta: %v", err)
	}
	return &w, nil
}

// ppmvToUgKg returns a multiplier to convert a concentration in
// ppmv dry air to a mass fraction [micrograms per kilogram dry air]
// for a chemical species with the given molecular weight in g/mol.
func ppmvToUgKg(mw float64) float64 {
	return mw * 1000.0 / MWa
}

func (w *WRFChem) read(varName string) NextData {
	return nextDataNCF(w.wrfOut, wrfFormat, varName, w.start, w.end, w.recordDelta, w.fileDelta, readNCF, w.msgChan)
}

func (w *WRFChem) readGroupAlt(varGroup map[string]float64) NextData {
	return nextDataGroupAltNCF(w.wrfOut, wrfFormat, varGroup, w.ALT(), w.start, w.end, w.recordDelta, w.fileDelta, readNCF, w.msgChan)
}

func (w *WRFChem) readGroup(varGroup map[string]float64) NextData {
	return nextDataGroupNCF(w.wrfOut, wrfFormat, varGroup, w.start, w.end, w.recordDelta, w.fileDelta, readNCF, w.msgChan)
}

// Nx helps fulfill the Preprocessor interface by returning
// the number of grid cells in the West-East direction.
func (w *WRFChem) Nx() (int, error) {
	f, ff, err := ncfFromTemplate(w.wrfOut, wrfFormat, w.start)
	if err != nil {
		return -1, fmt.Errorf("nx: %v", err)
	}
	defer f.Close()
	return ff.Header.Lengths("ALT")[3], nil
}

// Ny helps fulfill the Preprocessor interface by returning
// the number of grid cells in the South-North direction.
func (w *WRFChem) Ny() (int, error) {
	f, ff, err := ncfFromTemplate(w.wrfOut, wrfFormat, w.start)
	if err != nil {
		return -1, fmt.Errorf("ny: %v", err)
	}
	defer f.Close()
	return ff.Header.Lengths("ALT")[2], nil
}

// Nz helps fulfill the Preprocessor interface by returning
// the number of grid cells in the below-above direction.
func (w *WRFChem) Nz() (int, error) {
	f, ff, err := ncfFromTemplate(w.wrfOut, wrfFormat, w.start)
	if err != nil {
		return -1, fmt.Errorf("nz: %v", err)
	}
	defer f.Close()
	return ff.Header.Lengths("ALT")[1], nil
}

// PBLH helps fulfill the Preprocessor interface by returning
// planetary boundary layer height [m].
func (w *WRFChem) PBLH() NextData { return w.read("PBLH") }

// Height helps fulfill the Preprocessor interface by returning
// layer heights above ground level calculated based on geopotential height.
// For more information, refer to
// http://www.openwfm.org/wiki/How_to_interpret_WRF_variables.
func (w *WRFChem) Height() NextData {
	// ph is perturbation geopotential height [m2/s].
	phFunc := w.read("PH")
	// phb is baseline geopotential height [m2/s].
	phbFunc := w.read("PHB")
	return func() (*sparse.DenseArray, error) {
		ph, err := phFunc()
		if err != nil {
			return nil, err
		}
		phb, err := phbFunc()
		if err != nil {
			return nil, err
		}
		return geopotentialToHeight(ph, phb), nil
	}
}

func geopotentialToHeight(ph, phb *sparse.DenseArray) *sparse.DenseArray {
	layerHeights := sparse.ZerosDense(ph.Shape...)
	for k := 0; k < ph.Shape[0]; k++ {
		for j := 0; j < ph.Shape[1]; j++ {
			for i := 0; i < ph.Shape[2]; i++ {
				h := (ph.Get(k, j, i) + phb.Get(k, j, i) -
					ph.Get(0, j, i) - phb.Get(0, j, i)) / g // m
				layerHeights.Set(h, k, j, i)
			}
		}
	}
	return layerHeights
}

// ALT helps fulfill the Preprocessor interface by returning
// inverse air density [m3/kg].
func (w *WRFChem) ALT() NextData { return w.read("ALT") }

// U helps fulfill the Preprocessor interface by returning
// West-East wind speed [m/s].
func (w *WRFChem) U() NextData { return w.read("U") }

// V helps fulfill the Preprocessor interface by returning
// South-North wind speed [m/s].
func (w *WRFChem) V() NextData { return w.read("V") }

// W helps fulfill the Preprocessor interface by returning
// below-above wind speed [m/s].
func (w *WRFChem) W() NextData { return w.read("W") }

// AVOC helps fulfill the Preprocessor interface.
func (w *WRFChem) AVOC() NextData { return w.readGroupAlt(w.aVOC) }

// BVOC helps fulfill the Preprocessor interface.
func (w *WRFChem) BVOC() NextData { return w.readGroupAlt(w.bVOC) }

// NOx helps fulfill the Preprocessor interface.
func (w *WRFChem) NOx() NextData { return w.readGroupAlt(w.nox) }

// SOx helps fulfill the Preprocessor interface.
func (w *WRFChem) SOx() NextData { return w.readGroupAlt(w.sox) }

// NH3 helps fulfill the Preprocessor interface.
func (w *WRFChem) NH3() NextData { return w.readGroupAlt(w.nh3) }

// ASOA helps fulfill the Preprocessor interface.
func (w *WRFChem) ASOA() NextData { return w.readGroup(w.aSOA) }

// BSOA helps fulfill the Preprocessor interface.
func (w *WRFChem) BSOA() NextData { return w.readGroup(w.bSOA) }

// PNO helps fulfill the Preprocessor interface.
func (w *WRFChem) PNO() NextData { return w.readGroup(w.pNO) }

// PS helps fulfill the Preprocessor interface.
func (w *WRFChem) PS() NextData { return w.readGroup(w.pS) }

// PNH helps fulfill the Preprocessor interface.
func (w *WRFChem) PNH() NextData { return w.readGroup(w.pNH) }

// TotalPM25 helps fulfill the Preprocessor interface.
func (w *WRFChem) TotalPM25() NextData { return w.readGroup(w.totalPM25) }

// SurfaceHeatFlux helps fulfill the Preprocessor interface
// by returning heat flux at the surface [W/m2].
func (w *WRFChem) SurfaceHeatFlux() NextData { return w.read("HFX") }

// UStar helps fulfill the Preprocessor interface
// by returning friction velocity [m/s].
func (w *WRFChem) UStar() NextData { return w.read("UST") }

// T helps fulfill the Preprocessor interface by
// returning temperature [K].
func (w *WRFChem) T() NextData {
	thetaFunc := w.read("T") // perturbation potential temperature [K]
	pFunc := w.P()           // Pressure [Pa]
	return wrfTemperatureConvert(thetaFunc, pFunc)
}

func wrfTemperatureConvert(thetaFunc, pFunc NextData) NextData {
	return func() (*sparse.DenseArray, error) {
		thetaPerturb, err := thetaFunc() // perturbation potential temperature [K]
		if err != nil {
			return nil, err
		}
		p, err := pFunc() // Pressure [Pa]
		if err != nil {
			return nil, err
		}

		T := sparse.ZerosDense(thetaPerturb.Shape...)
		for i, tp := range thetaPerturb.Elements {
			T.Elements[i] = thetaPerturbToTemperature(tp, p.Elements[i])
		}
		return T, nil
	}
}

// thetaPerturbToTemperature converts perburbation potential temperature
// to ambient temperature for the given pressure (p [Pa]).
func thetaPerturbToTemperature(thetaPerturb, p float64) float64 {
	const (
		po    = 101300. // Pa, reference pressure
		kappa = 0.2854  // related to von karman's constant
	)
	pressureCorrection := math.Pow(p/po, kappa)
	// potential temperature, K
	θ := thetaPerturb + 300.
	// Ambient temperature, K
	return θ * pressureCorrection
}

// P helps fulfill the Preprocessor interface
// by returning pressure [Pa].
func (w *WRFChem) P() NextData {
	pbFunc := w.read("PB") // baseline pressure [Pa]
	pFunc := w.read("P")   // perturbation pressure [Pa]
	return wrfPressureConvert(pFunc, pbFunc)
}

func wrfPressureConvert(pFunc, pbFunc NextData) NextData {
	return func() (*sparse.DenseArray, error) {
		pb, err := pbFunc() // baseline pressure [Pa]
		if err != nil {
			return nil, err
		}
		p, err := pFunc() // perturbation pressure [Pa]
		if err != nil {
			return nil, err
		}
		P := pb.Copy()
		P.AddDense(p)
		return P, nil
	}
}

// HO helps fulfill the Preprocessor interface
// by returning hydroxyl radical concentration [ppmv].
func (w *WRFChem) HO() NextData { return w.read("OH") }

// H2O2 helps fulfill the Preprocessor interface
// by returning hydrogen peroxide concentration [ppmv].
func (w *WRFChem) H2O2() NextData { return w.read("HO2H") }

// SeinfeldLandUse helps fulfill the Preprocessor interface
// by returning land use categories as
// specified in github.com/ctessum/atmos/seinfeld.
func (w *WRFChem) SeinfeldLandUse() NextData {
	luFunc := w.read("LU_INDEX") // USGS land use index
	return wrfSeinfeldLandUse(luFunc)
}

func wrfSeinfeldLandUse(luFunc NextData) NextData {
	return func() (*sparse.DenseArray, error) {
		lu, err := luFunc() // USGS land use index
		if err != nil {
			return nil, err
		}
		o := sparse.ZerosDense(lu.Shape...)
		for j := 0; j < lu.Shape[0]; j++ {
			for i := 0; i < lu.Shape[1]; i++ {
				o.Set(float64(USGSseinfeld[f2i(lu.Get(j, i))-1]), j, i)
			}
		}
		return o, nil
	}
}

// USGSseinfeld lookup table to go from USGS land classes to land classes for
// particle dry deposition.
var USGSseinfeld = []seinfeld.LandUseCategory{
	seinfeld.Evergreen, //'Evergreen Needleleaf Forest'
	seinfeld.Deciduous, //'Evergreen Broadleaf Forest'
	seinfeld.Evergreen, //'Deciduous Needleleaf Forest'
	seinfeld.Deciduous, //'Deciduous Broadleaf Forest'
	seinfeld.Deciduous, //'Mixed Forest'
	seinfeld.Shrubs,    //'Closed Shrubland'
	seinfeld.Shrubs,    //'Open Shrubland'
	seinfeld.Shrubs,    //'Woody Savanna'
	seinfeld.Grass,     //'Savanna'
	seinfeld.Grass,     //'Grassland'
	seinfeld.Grass,     //'Permanent Wetland'
	seinfeld.Grass,     //'Cropland'
	seinfeld.Desert,    //'Urban and Built-Up'
	seinfeld.Grass,     //'Cropland / Natural Veg. Mosaic'
	seinfeld.Desert,    //'Permanent Snow'
	seinfeld.Desert,    //'Barren / Sparsely Vegetated'
	seinfeld.Desert,    //'IGBP Water'
	seinfeld.Desert,    //'Unclassified'
	seinfeld.Desert,    //'Fill Value'
	seinfeld.Desert,    //'Unclassified'
	seinfeld.Desert,    //'Open Water'
	seinfeld.Desert,    //'Perennial Ice/Snow'
	seinfeld.Desert,    //'Developed Open Space'
	seinfeld.Desert,    //'Developed Low Intensity'
	seinfeld.Desert,    //'Developed Medium Intensity'
	seinfeld.Desert,    //'Developed High Intensity'
	seinfeld.Desert,    //'Barren Land'
	seinfeld.Deciduous, //'Deciduous Forest'
	seinfeld.Evergreen, //'Evergreen Forest'
	seinfeld.Deciduous, //'Mixed Forest'
	seinfeld.Shrubs,    //'Dwarf Scrub'
	seinfeld.Shrubs,    //'Shrub/Scrub'
	seinfeld.Grass,     //'Grassland/Herbaceous'
	seinfeld.Grass,     //'Sedge/Herbaceous'
	seinfeld.Desert,    //'Lichens'
	seinfeld.Desert,    //'Moss'
	seinfeld.Grass,     //'Pasture/Hay'
	seinfeld.Grass,     //'Cultivated Crops'
	seinfeld.Deciduous, //'Woody Wetland'
	seinfeld.Grass,     //'Emergent Herbaceous Wetland'
}

// WeselyLandUse helps fulfill the Preprocessor interface
// by returning land use categories as
// specified in github.com/ctessum/atmos/wesely1989.
func (w *WRFChem) WeselyLandUse() NextData {
	luFunc := w.read("LU_INDEX") // USGS land use index
	return wrfWeselyLandUse(luFunc)
}

func wrfWeselyLandUse(luFunc NextData) NextData {
	return func() (*sparse.DenseArray, error) {
		lu, err := luFunc() // USGS land use index
		if err != nil {
			return nil, err
		}
		o := sparse.ZerosDense(lu.Shape...)
		for j := 0; j < lu.Shape[0]; j++ {
			for i := 0; i < lu.Shape[1]; i++ {
				o.Set(float64(USGSwesely[f2i(lu.Get(j, i))-1]), j, i)
			}
		}
		return o, nil
	}
}

// USGSwesely lookup table to go from USGS land classes to land classes for
// gas dry deposition.
var USGSwesely = []wesely1989.LandUseCategory{
	wesely1989.Coniferous,  //'Evergreen Needleleaf Forest'
	wesely1989.Deciduous,   //'Evergreen Broadleaf Forest'
	wesely1989.Coniferous,  //'Deciduous Needleleaf Forest'
	wesely1989.Deciduous,   //'Deciduous Broadleaf Forest'
	wesely1989.MixedForest, //'Mixed Forest'
	wesely1989.RockyShrubs, //'Closed Shrubland'
	wesely1989.RockyShrubs, //'Open Shrubland'
	wesely1989.RockyShrubs, //'Woody Savanna'
	wesely1989.Range,       //'Savanna'
	wesely1989.Range,       //'Grassland'
	wesely1989.Wetland,     //'Permanent Wetland'
	wesely1989.RangeAg,     //'Cropland'
	wesely1989.Urban,       //'Urban and Built-Up'
	wesely1989.RangeAg,     //'Cropland / Natural Veg. Mosaic'
	wesely1989.Barren,      //'Permanent Snow'
	wesely1989.Barren,      //'Barren / Sparsely Vegetated'
	wesely1989.Water,       //'IGBP Water'
	wesely1989.Barren,      //'Unclassified'
	wesely1989.Barren,      //'Fill Value'
	wesely1989.Barren,      //'Unclassified'
	wesely1989.Water,       //'Open Water'
	wesely1989.Barren,      //'Perennial Ice/Snow'
	wesely1989.Urban,       //'Developed Open Space'
	wesely1989.Urban,       //'Developed Low Intensity'
	wesely1989.Urban,       //'Developed Medium Intensity'
	wesely1989.Urban,       //'Developed High Intensity'
	wesely1989.Barren,      //'Barren Land'
	wesely1989.Deciduous,   //'Deciduous Forest'
	wesely1989.Coniferous,  //'Evergreen Forest'
	wesely1989.MixedForest, //'Mixed Forest'
	wesely1989.RockyShrubs, //'Dwarf Scrub'
	wesely1989.RockyShrubs, //'Shrub/Scrub'
	wesely1989.Range,       //'Grassland/Herbaceous'
	wesely1989.Range,       //'Sedge/Herbaceous'
	wesely1989.Barren,      //'Lichens'
	wesely1989.Barren,      //'Moss'
	wesely1989.RangeAg,     //'Pasture/Hay'
	wesely1989.RangeAg,     //'Cultivated Crops'
	wesely1989.Wetland,     //'Woody Wetland'
	wesely1989.Wetland,     //'Emergent Herbaceous Wetland'
}

// Z0 helps fulfill the Preprocessor interface by
// returning roughness length.
func (w *WRFChem) Z0() NextData {
	LUIndexFunc := w.read("LU_INDEX") //USGS land use index
	return wrfZ0(LUIndexFunc)
}

func wrfZ0(LUIndexFunc NextData) NextData {
	return func() (*sparse.DenseArray, error) {
		luIndex, err := LUIndexFunc()
		if err != nil {
			return nil, err
		}
		zo := sparse.ZerosDense(luIndex.Shape...)
		for i, lu := range luIndex.Elements {
			zo.Elements[i] = USGSz0[f2i(lu)-1] // roughness length [m]
		}
		return zo, nil
	}
}

// USGSz0 holds Roughness lengths for USGS land classes ([m]), from WRF file
// VEGPARM.TBL.
var USGSz0 = []float64{.50, .50, .50, .50, .35, .03, .035, .03, .15, .11,
	.30, .10, .50, .095, .001, .01, .0001, 999., 999., 999.,
	.0001, .001, .50, .70, 1.5, 2.0, .01, .50, .50, .35,
	.025, .03, .11, .20, .01, .01, .10, .06, .40, .20}

// QRain helps fulfill the Preprocessor interface by
// returning rain mass fraction.
func (w *WRFChem) QRain() NextData { return w.read("QRAIN") }

// CloudFrac helps fulfill the Preprocessor interface
// by returning the fraction of each grid cell filled
// with clouds [volume/volume].
func (w *WRFChem) CloudFrac() NextData { return w.read("CLDFRA") }

// QCloud helps fulfill the Preprocessor interface by returning
// the mass fraction of cloud water in each grid cell [mass/mass].
func (w *WRFChem) QCloud() NextData { return w.read("QCLOUD") }

// RadiationDown helps fulfill the Preprocessor interface by returning
// total downwelling radiation at ground level [W/m2].
func (w *WRFChem) RadiationDown() NextData {
	swDownFunc := w.read("SWDOWN") // downwelling short wave radiation at ground level [W/m2]
	glwFunc := w.read("GLW")       // downwelling long wave radiation at ground level [W/m2]
	return wrfRadiationDown(swDownFunc, glwFunc)
}

func wrfRadiationDown(swDownFunc, glwFunc NextData) NextData {
	return func() (*sparse.DenseArray, error) {
		swDown, err := swDownFunc() // downwelling short wave radiation at ground level [W/m2]
		if err != nil {
			return nil, err
		}
		glw, err := glwFunc() // downwelling long wave radiation at ground level [W/m2]
		if err != nil {
			return nil, err
		}
		rad := swDown.Copy()
		rad.AddDense(glw)
		return rad, nil
	}
}

// SWDown helps fulfill the Preprocessor interface by returning
// downwelling short wave radiation at ground level [W/m2].
func (w *WRFChem) SWDown() NextData { return w.read("SWDOWN") }

// GLW helps fulfill the Preprocessor interface by returning
// downwelling long wave radiation at ground level [W/m2].
func (w *WRFChem) GLW() NextData { return w.read("GLW") }
