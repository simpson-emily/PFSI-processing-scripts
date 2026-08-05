
// Supplementary information to: Simpson et al. 2026: Fire-Moisture Interactions Shape Recovery Trajectories Across Vegetation Communities, Global Change Biology
// ** SCRIPT 1 **

// -----------------------------------------------------------------------------
// 1. VEGETATION REMAPPING 
// vegRaw = NSW State Type Vegetation Map (DCCEEW): https://datasets.seed.nsw.gov.au/dataset/nsw-state-vegetation-type-map
// aoi = Greater Blue Mountains World Heritage Area
// vegTable = set of raw pixel values and their corresponding community codes
// -----------------------------------------------------------------------------

// Community 10 = Not Classified
var vegTableFiltered = vegTable.filter(ee.Filter.neq('community', 10));

// Create unique community codes
var uniqueCommunities = vegTableFiltered.aggregate_array('community').distinct();
var commCodes = ee.List.sequence(0, uniqueCommunities.length().subtract(1));
var commDict = ee.Dictionary.fromLists(uniqueCommunities, commCodes);

var vegTableWithCodes = vegTableFiltered.map(function(f) {
  var original = f.get('community');
  var code = commDict.get(original);
  return f.set('comm_code', code);
});

// Mask the raw layer.
var vegRawMasked = vegRaw.updateMask(vegRaw.neq(10));

// Remap raw_value to the community code
var fromVals = vegTableWithCodes.aggregate_array('raw_value');
var toVals   = vegTableWithCodes.aggregate_array('comm_code');
var vegLayer = vegRawMasked.remap(fromVals, toVals)
                     .toInt()
                     .rename('VegType');

// Merge Dry Sclerophyll Forest and Wet Sclerophyll Forest types
vegLayer = vegLayer
  .where(vegLayer.eq(3).or(vegLayer.eq(4)),  34)
  .where(vegLayer.eq(15).or(vegLayer.eq(16)), 1516)
  .rename('VegType');

vegLayer = vegLayer.updateMask(vegLayer.neq(10));

// Update vegLayer once it is fully built and masked
vegLayer = vegLayer.updateMask(vegLayer.neq(10));

// Export to use when sampling
Export.image.toDrive({
   image:       vegLayer,
   description: 'VegLayer',
   folder:      'your_folder_here',       
   fileNamePrefix: 'VegLayer', 
   region:      aoi,                 
   scale:       30,                  
   maxPixels:   1e9,                 
   crs:         'EPSG:4326'          
  });


// -----------------------------------------------------------------------------
// 2. LOAD DATA, CLOUD MASKING AND SCALING
// Process a season at a time. Comment out unused block depending on what season is being processed.
// -----------------------------------------------------------------------------

// For Spring, Winter, Autumn access data as:

var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterDate('2014-03-01', '2024-5-31')
  .filter(ee.Filter.calendarRange(3, 5, 'month'))
  .filterBounds(aoi);
  
// For Summer, we need Dec of the previous year:

var l8sri = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(aoi)
  .filterDate('2013-12-01', '2024-02-29');  

// Define the years for which you want to extract summer seasons
var years = ee.List.sequence(2014, 2024);

// Function extracts Dec (previous year), Jan, Feb (current year) images
var getSummerImages = function(year) {
  year = ee.Number(year);
  var decPrevYear = l8sri.filterDate(ee.Date.fromYMD(year.subtract(1), 12, 1), ee.Date.fromYMD(year, 1, 1));
  var janFebThisYear = l8sri.filterDate(ee.Date.fromYMD(year, 1, 1), ee.Date.fromYMD(year, 3, 1));
  return decPrevYear.merge(janFebThisYear)
    .map(function(img) {
      return img.set('summer_year', year);
    });
};

var l8 = ee.ImageCollection(
  years.iterate(function(year, acc) {
    var collection = getSummerImages(year);
    return ee.ImageCollection(acc).merge(collection);
  }, ee.ImageCollection([]))
);

// Correct for clouds and scale bands.
function maskAndScale(img) {
  var qa = img.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1<<3).eq(0)
           .and(qa.bitwiseAnd(1<<4).eq(0));
  var scaled = img.select(
      ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'],
      ['Blue','Green','Red','NIR','SWIR1','SWIR2']
    )
    .multiply(0.0000275).add(-0.2);
  return img.addBands(scaled, null, true)
            .updateMask(mask)
            .copyProperties(img, ['system:time_start']);
}

var prepped = l8.map(maskAndScale);

// -----------------------------------------------------------------------------
// 3. LOAD TERRAIN MODELS. CALCULATE CANOPY HEIGHT MODEL.
// -----------------------------------------------------------------------------

var dem     = ee.Image('USGS/SRTMGL1_003').clip(aoi); // ground elevation map
var terrain = ee.Terrain.products(dem); // GEE's inbuilt terrain algorithm
var slope   = terrain.select('slope'); // calculates how steep each pixel is in degrees
var aspect  = terrain.select('aspect'); // calculates compass direction of each pixel

var aw3dTiles = ee.ImageCollection('JAXA/ALOS/AW3D30/V4_1') 
  .filterBounds(aoi)             
  .select('DSM'); 

var dsm30 = aw3dTiles
  .mosaic()                      
  .rename('DSM30')
  .clip(aoi);
var dtm30 = dem;  

var chm30 = dsm30.subtract(dtm30) // computes Canopy Height Model
  .max(0) // removes negs ('below ground' = invalid heights)
  .rename('CHM')
  .clip(aoi);
  
// -----------------------------------------------------------------------------
// 4. SOLAR & TERRAIN GEOMETRY
// -----------------------------------------------------------------------------

function addGeometry(img) {
  var sunEl = ee.Number(img.get('SUN_ELEVATION')); // suns elevation above horizon
  var sunAz = ee.Number(img.get('SUN_AZIMUTH')); // compass direction of sun

  var ziNum    = ee.Number(90).subtract(sunEl).multiply(Math.PI/180);
  var relAzNum = sunAz.multiply(-1).multiply(Math.PI/180); //the angle between the sun and sensor direction (calc by x -1)

  var ziImg    = ee.Image.constant(ziNum).toFloat().rename('zi');
  var relAzImg = ee.Image.constant(relAzNum).toFloat().rename('relAz');

  var sRad = slope.multiply(Math.PI/180).toFloat().rename('slopeRad');
  var aRad = aspect.multiply(Math.PI/180).toFloat().rename('aspectRad');

  return img.addBands([ziImg, relAzImg, sRad, aRad]);
}

// -----------------------------------------------------------------------------
// 5. CROWN RADIUS (R) & HEIGHT‐TO‐RADIUS (H)
// -----------------------------------------------------------------------------

var vegType = vegLayer.clip(aoi);

var Rmap = vegType.remap( 
  [1516, 34, 5, 6, 8, 9, 11], // vegetation community codes
  [  10,  8, 4, 1, 3, 2, 12], // radius: values are averages from combination of in-field observation and literature
  10.0)
  .rename('R'); 

var Hmap = chm30.divide(Rmap) 
  .clamp(0.05, 2.0)
  .rename('H'); // H = how tall each canopy is compared to its width
  
// To note:
// 1516 = Wet Sclerophyll Forest
// 34 = Dry Sclerophyll Forest
// 5 = Forested Wetland
// 6 = Freshwater Wetland
// 8 = Grassy Woodland
// 9 = Heathland
// 11 = Rainforest
  
// -----------------------------------------------------------------------------
// 6. CALCULATE BRDF KERNELS (KVOL / KGEO)
// -----------------------------------------------------------------------------

// Calculate Kvol (volume-scattering kernel)

function rossThick(zi, zv, relAz) { 
  var cosZi = zi.cos(), cosZv = zv.cos();
  var term1 = zi.subtract(zv).abs()
               .multiply(-1)
               .add(Math.PI/2).cos();
  return term1
    .subtract(cosZi.multiply(cosZv))
    .divide(cosZi.multiply(cosZv))
    .rename('Kvol'); 
}

// Calculate Kgeo. How much of the canopy is shadowed versus lit for any sun/view geometry?

function liSparse(zi, zv, relAz, slopeRad, aspectRad, H) { // tilt sun and view angles based on canopy height (taller canopies change path of light)
  var th_i_p = zi.tan().multiply(H).atan(); // slope of incoming sunlight
  var th_v_p = zv.tan().multiply(H).atan();
  var ci_p = th_i_p.cos(), cv_p = th_v_p.cos();
  var si_p = th_i_p.sin(), sv_p = th_v_p.sin();
  var cos_t = relAz.cos()
      .multiply(si_p).multiply(sv_p)
      .add(ci_p.multiply(cv_p));
  var t = cos_t.acos();
  
  var sec_i = ci_p.pow(-1), sec_v = cv_p.pow(-1); //how much of the canopy 'shadows' itself
  var O = t.subtract(t.sin().multiply(t.cos()))
           .divide(Math.PI)
           .multiply(sec_i.add(sec_v));
  var cz_i = zi.cos(), cz_v = zv.cos();
  
  var Kgeo = O.subtract( // build Kgeo (geometric kernel)
               cz_i.multiply(cz_v).divide(cz_i.add(cz_v))
             ).rename('Kgeo');
  return Kgeo;
}

function addKernels(img) {
  var zi    = img.select('zi'),
      zv    = ee.Image.constant(0),
      relAz = img.select('relAz'),
      slp   = img.select('slopeRad'),
      asp   = img.select('aspectRad'),
      H     = img.select('H');
  var Kvol = rossThick(zi, zv, relAz),
      Kgeo = liSparse(zi, zv, relAz, slp, asp, H);
  return img.addBands([Kvol, Kgeo]);
}

// -----------------------------------------------------------------------------
// 7. BRDF + TOPOGRAPHIC CORRECTION 
// -----------------------------------------------------------------------------

// Constants from: Roy et al. 2016 "A general method to normalize Landsat reflectance data to nadir BRDF adjusted reflectance"

var F_iso  = ee.Image.constant([0.0774,0.1306,0.1690,0.3146,0.1884,0.1151]); 
var F_vol  = ee.Image.constant([0.0372,0.0585,0.0762,0.1589,0.0679,0.0343]); 
var F_geo  = ee.Image.constant([0.0079,0.0174,0.0222,0.0388,0.0160,0.0086]); 
var cosRef = ee.Image.constant(Math.cos(45 * Math.PI/180)); 

function applyBRDFTopo(img) {
  var refl     = img.select(['Blue','Green','Red','NIR','SWIR1','SWIR2']);
  var Kvol     = img.select('Kvol');
  var Kgeo     = img.select('Kgeo');
  var zi       = img.select('zi');
  var slp      = img.select('slopeRad');
  var relAz    = img.select('relAz');
  
  var cosTheta = zi.cos().multiply(slp.cos())
      .add(zi.sin().multiply(slp.sin()).multiply(relAz.cos()));
  var brdfObs  = F_iso
                 .add(F_vol.multiply(Kvol))
                 .add(F_geo.multiply(Kgeo));
  var cFac     = F_iso.divide(brdfObs);
  var nbar     = refl.multiply(cFac)
                    .multiply(cosRef.divide(cosTheta))
                    .rename([
                      'NBAR_Blue','NBAR_Green','NBAR_Red',
                      'NBAR_NIR','NBAR_SWIR1','NBAR_SWIR2'
                    ]);
  return img.addBands(nbar, null, true);
}

// Build the NBAR‐corrected collection
var correctedCollection = prepped
  .map(addGeometry)
  .map(function(img) {
    return img.addBands([chm30.rename('CHM'), Rmap, Hmap]);
  })
  .map(addKernels)
  .map(applyBRDFTopo);

// -----------------------------------------------------------------------------
// 8. COMPOSITE IMAGES
// -----------------------------------------------------------------------------

// Build yearly composites for season of interest 
var years = ee.List.sequence(2014, 2024);
var yearlyComposites = ee.ImageCollection(
  years.map(function(y) {
    return correctedCollection
      .filter(ee.Filter.calendarRange(y, y, 'year'))
      .sort('system:time_start')
      .median()
      .clip(aoi)
      .set('year', y);
  })
);

var bestImage2014 = yearlyComposites
  .filter(ee.Filter.eq('year', 2014))
  .first();
var bestImage2015 = yearlyComposites
  .filter(ee.Filter.eq('year', 2015))
  .first();
var bestImage2016 = yearlyComposites
  .filter(ee.Filter.eq('year', 2016))
  .first();
var bestImage2017 = yearlyComposites
  .filter(ee.Filter.eq('year', 2017))
  .first();
var bestImage2018 = yearlyComposites
  .filter(ee.Filter.eq('year', 2018))
  .first();
var bestImage2019 = yearlyComposites
  .filter(ee.Filter.eq('year', 2019))
  .first();
var bestImage2020 = yearlyComposites
  .filter(ee.Filter.eq('year', 2020))
  .first();
var bestImage2021 = yearlyComposites
  .filter(ee.Filter.eq('year', 2021))
  .first();
var bestImage2022 = yearlyComposites
  .filter(ee.Filter.eq('year', 2022))
  .first();
var bestImage2023 = yearlyComposites
  .filter(ee.Filter.eq('year', 2023))
  .first();
var bestImage2024 = yearlyComposites
  .filter(ee.Filter.eq('year', 2024))
  .first();

// -----------------------------------------------------------------------------
// 9. CALCULATE NBR2
// -----------------------------------------------------------------------------

function calculateNBR2(image) {
  var swir1       = image.select('NBAR_SWIR1');
  var swir2       = image.select('NBAR_SWIR2');
  var denominator = swir1.add(swir2);

  var nbr2_fraction = swir1
    .subtract(swir2)
    .divide(denominator)
    .rename('NBR2_frac');

  // Scale
  var nbr2_int16 = nbr2_fraction
    .multiply(10000)       
    .round()               
    .toShort()             
    .max(-10000)           
    .min(10000)            
    .rename('NBR2');       

  return image
    .addBands(nbr2_int16, null, true)
    .copyProperties(image, ['system:time_start']);
}

var NBR2014 = calculateNBR2(bestImage2014);
var NBR2015 = calculateNBR2(bestImage2015);
var NBR2016 = calculateNBR2(bestImage2016);
var NBR2017 = calculateNBR2(bestImage2017);
var NBR2018 = calculateNBR2(bestImage2018);
var NBR2019 = calculateNBR2(bestImage2019);
var NBR2020 = calculateNBR2(bestImage2020);
var NBR2021 = calculateNBR2(bestImage2021);
var NBR2022 = calculateNBR2(bestImage2022);
var NBR2023 = calculateNBR2(bestImage2023);
var NBR2024 = calculateNBR2(bestImage2024);

// -----------------------------------------------------------------------------
// 10. CALCULATE POST-FIRE STABILITY INDEX
// Gibson et al. 2022, "The post-fire stability index; a new approach to monitoring 
// post-fire recovery by satellite imagery"
// -----------------------------------------------------------------------------

function calculatePFSI(current, previous) {
  current  = ee.Image(current).select('NBR2');
  previous = ee.Image(previous).select('NBR2');

  var pfs = current
    .subtract(previous)
    .divide(
      previous
        .divide(1000)
        .abs()
        .sqrt()
    )
    .rename('PFSI');

  return pfs.copyProperties(current, ['system:time_start']);
}

var PFSI2015 = calculatePFSI(NBR2015, NBR2014);
var PFSI2016 = calculatePFSI(NBR2016, NBR2015);
var PFSI2017 = calculatePFSI(NBR2017, NBR2016);
var PFSI2018 = calculatePFSI(NBR2018, NBR2017);
var PFSI2019 = calculatePFSI(NBR2019, NBR2018);
var PFSI2020 = calculatePFSI(NBR2020, NBR2019);
var PFSI2021 = calculatePFSI(NBR2021, NBR2020);
var PFSI2022 = calculatePFSI(NBR2022, NBR2021);
var PFSI2023 = calculatePFSI(NBR2023, NBR2022);
var PFSI2024 = calculatePFSI(NBR2024, NBR2023);


// OPTIONAL(to visualise): Can classify raw values based on Table 3 in Gibson et al. (2022)

function classifyPFSI(pfsiImage) {
  var classified = pfsiImage
    .where(pfsiImage.gt(500), 13)                            // extreme increase
    .where(pfsiImage.gt(400).and(pfsiImage.lte(500)), 12)    // very large increase
    .where(pfsiImage.gt(300).and(pfsiImage.lte(400)), 11)    // large increase
    .where(pfsiImage.gt(200).and(pfsiImage.lte(300)), 10)    // moderate increase
    .where(pfsiImage.gt(100).and(pfsiImage.lte(200)), 9)     // small increase
    .where(pfsiImage.gt(50).and(pfsiImage.lte(100)), 8)      // very small increase
    .where(pfsiImage.gte(-50).and(pfsiImage.lte(50)), 7)     // stable range
    .where(pfsiImage.gt(-100).and(pfsiImage.lt(-50)), 6)     // very small decrease
    .where(pfsiImage.gt(-200).and(pfsiImage.lte(-100)), 5)   // small decrease
    .where(pfsiImage.gt(-300).and(pfsiImage.lte(-200)), 4)   // moderate decrease
    .where(pfsiImage.gt(-400).and(pfsiImage.lte(-300)), 3)   // large decrease
    .where(pfsiImage.gt(-500).and(pfsiImage.lte(-400)), 2)   // very large decrease
    .where(pfsiImage.lt(-500), 1);                           // extreme decrease

  return classified;
}

var classifiedPFSI_2015 = classifyPFSI(ee.Image(PFSI2015).select('PFSI'));
var classifiedPFSI_2016 = classifyPFSI(ee.Image(PFSI2016).select('PFSI'));
var classifiedPFSI_2017 = classifyPFSI(ee.Image(PFSI2017).select('PFSI'));
var classifiedPFSI_2018 = classifyPFSI(ee.Image(PFSI2018).select('PFSI'));
var classifiedPFSI_2019 = classifyPFSI(ee.Image(PFSI2019).select('PFSI'));
var classifiedPFSI_2020 = classifyPFSI(ee.Image(PFSI2020).select('PFSI'));
var classifiedPFSI_2021 = classifyPFSI(ee.Image(PFSI2021).select('PFSI'));
var classifiedPFSI_2022 = classifyPFSI(ee.Image(PFSI2022).select('PFSI'));
var classifiedPFSI_2023 = classifyPFSI(ee.Image(PFSI2023).select('PFSI'));
var classifiedPFSI_2024 = classifyPFSI(ee.Image(PFSI2024).select('PFSI'));


// Choose colours
var classificationVis = {
  min: 1,
  max: 13,
  palette: ['#9e1b00', '#bf360c', '#e64a19','#ff5f2d', '#ff8c4e', '#ffb37a', '#f5f5dc', 
  '#b8e6ea', '#7cccd6', '#4ea3c3', '#2c89bc', '#1b6ca8','#004c6d'
]
};

// Add the classified layer to map to view

// Map.addLayer(classifiedPFSI_2015, classificationVis, 'Classified PFSI 2015');
// Map.addLayer(classifiedPFSI_2016, classificationVis, 'Classified PFSI 2016');
// Map.addLayer(classifiedPFSI_2017, classificationVis, 'Classified PFSI 2017');
// Map.addLayer(classifiedPFSI_2018, classificationVis, 'Classified PFSI 2018');
// Map.addLayer(classifiedPFSI_2019, classificationVis, 'Classified PFSI 2019');
// Map.addLayer(classifiedPFSI_2020, classificationVis, 'Classified PFSI 2020');
// Map.addLayer(classifiedPFSI_2021, classificationVis, 'Classified PFSI 2021');
// Map.addLayer(classifiedPFSI_2022, classificationVis, 'Classified PFSI 2022');
// Map.addLayer(classifiedPFSI_2023, classificationVis, 'Classified PFSI 2023');
// Map.addLayer(classifiedPFSI_2024, classificationVis, 'Classified PFSI 2024');


// To export images to google drive. One image exported for each season of each year:
// for example:

Export.image.toDrive({
  image: PFSI2016,
  description: '2016_Autumn',
  fileNamePrefix: '2016_Autumn',
  region: aoi,
  scale: 30, 
  crs: 'EPSG:4326', 
  maxPixels: 1e9,
  fileFormat: 'GeoTIFF'
});