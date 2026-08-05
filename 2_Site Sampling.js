
// Supplementary information to: Simpson et al. 2026: Fire-Moisture Interactions Shape Recovery Trajectories Across Vegetation Communities, Global Change Biology
// ** SCRIPT 2 **

// fireRaw = Fire Extent and Severity Map (native 10m, FESM)
// https://datasets.seed.nsw.gov.au/dataset/fire-extent-and-severity-mapping-fesm-2019-20
// vegLayerf = vegLayer from SCRIPT 1 (native 5m, SVTM)
// wha equivalent to aoi from SCRIPT 1

var vegLayers = vegLayerf.rename('VegType');
var vegLayer = vegLayers.updateMask(vegLayers.neq(0));

// -----------------------------------------------------------------------------
// 1. AGGREGATING VEG AND FIRE SEV
// -----------------------------------------------------------------------------

var targetProj = fireRaw.projection().atScale(30); 

// aggregate veg to 30m
var veg3x3 = vegLayer
  .reduceResolution({
    reducer: ee.Reducer.mode(), // chooses most common veg type
    maxPixels: 4096
  })
  .reproject({
    crs: targetProj
  })
  .rename('VegType_mode');

// aggregate fire sev to 30m
var fireLayer = fireRaw.round().toInt().rename("FireSeverity");

var fire3x3 = fireLayer
  .reduceResolution({
    reducer: ee.Reducer.mode(), // chooses most common fire sev
    maxPixels: 1024
  })
  .reproject({
    crs: targetProj
  })
  .rename('FireSeverity_mode');

// Combine layers for sampling
var combinedClass = veg3x3.multiply(10).add(fire3x3).toInt().rename('CombinedClass');

var baseStack = combinedClass
  .addBands(veg3x3)
  .addBands(fire3x3)
  .updateMask(fire3x3.neq(1));

// -----------------------------------------------------------------------------
// 2. CREATE/SAMPLE SITES (30 POINTS PER FIRE SEVERITY X VEG TYPE)
// -----------------------------------------------------------------------------

var combinedClassFreq = baseStack.select('CombinedClass').reduceRegion({
  reducer: ee.Reducer.frequencyHistogram(),
  geometry: wha,
  scale: 30,
  maxPixels: 1e13
}).get('CombinedClass');

var observedClasses = ee.Dictionary(combinedClassFreq).keys().map(ee.Number.parse);

var pointsPerClass = observedClasses.map(function(_) { return 30; });

// Generate candidate samples
var candidateSamples = baseStack.select('CombinedClass').stratifiedSample({
  numPoints: 0,
  classBand: 'CombinedClass',
  classValues: observedClasses,
  classPoints: pointsPerClass,
  region: wha,
  scale: 30,
  seed: 42,
  geometries: true,
  dropNulls: true
});

// Make sure no sampled points overlap. If they do, remove.

var bufferRadius = 45; 

// Assign dominant fire severity and veg type 

var bufferedSamples = candidateSamples.map(function(f) {
  var geom = f.geometry();

  // FireSeverity_mode from fire3x3 image
  var fsVal = fire3x3.reduceRegion({
    reducer:   ee.Reducer.first(),
    geometry:  geom,
    scale:     30,
    maxPixels: 1e9
  }).get('FireSeverity_mode');

  // VegType_mode from veg3x3 image
  var vegVal = veg3x3.reduceRegion({
    reducer:   ee.Reducer.first(),
    geometry:  geom,
    scale:     30,
    maxPixels: 1e9
  }).get('VegType_mode');

  return f.set({
    'FireSeverity_mode3x3': fsVal,
    'VegType_mode3x3'     : vegVal
  });
});

// Filter to keep non overlapping points
bufferedSamples = bufferedSamples
  .filter(ee.Filter.neq('FireSeverity_mode3x3', null))
  .filter(ee.Filter.neq('VegType_mode3x3',     null));

var distinctSamples = ee.List([]);
var sampleList = bufferedSamples.toList(bufferedSamples.size());

var filteredList = ee.List(sampleList.iterate(function(f, acc) {
  f = ee.Feature(f);
  acc = ee.List(acc);

  var fGeom = f.geometry().buffer(bufferRadius * 2);
  var intersects = ee.FeatureCollection(acc).filterBounds(fGeom).size();

  return ee.Algorithms.If(
    intersects.eq(0),
    acc.add(f),
    acc
  );
}, distinctSamples));

var masterSamples = ee.FeatureCollection(filteredList);

// Assign unique sample (i.e. site) IDs to all remaining points 

var size = masterSamples.size();
var fcList = masterSamples.toList(size);

var masterSamplesWithID = ee.FeatureCollection(
  ee.List.sequence(0, size.subtract(1)).map(function(i) {
    var f = ee.Feature(fcList.get(i));
    return f.set('sample_id', i);
  })
);

masterSamplesWithID = ee.FeatureCollection(masterSamplesWithID.toList(masterSamplesWithID.size()));

// -----------------------------------------------------------------------------
// 3. LOAD PFSI AND SM FILES. 
// Paths below are placeholders referencing PFSI composites created in SCRIPT 1, and 
// uploaded to authors GEE assets. Paths also reference SM composites generated 
// from rasters available publicly at: https://shiny.esoil.io/SMIPS/.
// -----------------------------------------------------------------------------

var pfsiImages = {
  '2015_Summer': ee.Image('your/path/here/2015_Summer_Annual'),
  '2015_Autumn': ee.Image('your/path/here/2015_Autumn_Annual'),
  '2015_Winter': ee.Image('your/path/here/2015_Winter_Annual'),
  '2015_Spring': ee.Image('your/path/here/2015_Spring_Annual'),

//.......

  '2024_Summer': ee.Image('your/path/here/2024_Summer_Annual'),
  '2024_Autumn': ee.Image('your/path/here/2024_Autumn_Annual'),
  '2024_Winter': ee.Image('your/path/here/2024_Winter_Annual'),
  '2024_Spring': ee.Image('your/path/here/2024_Spring_Annual')
};

var smImages = {
  '2014_Summer': ee.Image('your/path/here/2014_Summer_SM'),
  '2014_Autumn': ee.Image('your/path/here/2014_Autumn_SM'),
  '2014_Winter': ee.Image('your/path/here/2014_Winter_SM'),
  '2014_Spring': ee.Image('your/path/here/2014_Spring_SM'),
  
//.......

  '2024_Summer': ee.Image('your/path/here/2024_Summer_SM'),
  '2024_Autumn': ee.Image('your/path/here/2024_Autumn_SM'),
  '2024_Winter': ee.Image('your/path/here/2024_Winter_SM'),
  '2024_Spring': ee.Image('your/path/here/2024_Spring_SM')
};

// -----------------------------------------------------------------------------
// 4. SAMPLE PFSI AND SM AT SITES
// -----------------------------------------------------------------------------

// PFSI
var pfsi3x3Dict = ee.Dictionary(pfsiImages).map(function(key, img) {
  return ee.Image(img).rename([key]); 
});

// SM 
var sm3x3Dict = ee.Dictionary(smImages).map(function(key, img) {
  return ee.Image(img).rename([key]);
});

// Select ONLY those sample points which have valid PFSI & Soil Moisture for all year / season combinations

// Build list of valid year-season keys (i.e. create all possible combinations of year/season)
var pfsiKeys = Object.keys(pfsiImages);
var validKeys = pfsiKeys.filter(function(key) {
  return smImages.hasOwnProperty(key);
});

validKeys = ee.List(validKeys);

var validYearSeasonPairs = validKeys.map(function(key) {
  var parts = ee.String(key).split('_');
  return ee.List([ee.Number.parse(parts.get(0)), parts.get(1)]);
});


// Then check if points have valid values for PFSI and SM data for all given year/season keys
function hasValidDataForYS(point, year, season) {
  var key = ee.String(year).cat('_').cat(season);
  var pfsi = ee.Image(pfsi3x3Dict.get(key));
  var sm   = ee.Image(sm3x3Dict.get(key));

  var pfsiVal = pfsi.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: point.geometry(),
    scale: 30,
    maxPixels: 1e9
  }).values().get(0);

  var smVal = sm.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: point.geometry(),
    scale: 30,
    maxPixels: 1e9
  }).values().get(0);

  return ee.Algorithms.If(pfsiVal, ee.Algorithms.If(smVal, true, false), false);
}

function hasValidDataForAllYearSeasons(point) {
  return validYearSeasonPairs.iterate(function(pair, validSoFar) {
    var year = ee.Number(ee.List(pair).get(0));
    var season = ee.String(ee.List(pair).get(1));
    
    return ee.Algorithms.If(
      validSoFar,
      hasValidDataForYS(point, year, season),
      false
    );
  }, true);
}

// Creates column 'valid' in dataset. Values of 0 = not valid, 1 = valid, can keep
var flaggedSamples = masterSamplesWithID.map(function(f) {
  var isValid = hasValidDataForAllYearSeasons(f);
  return ee.Feature(f).set('valid', ee.Algorithms.If(isValid, 1, 0));
});

// Filter for valid points only
var filteredMasterSamples = flaggedSamples.filter(ee.Filter.eq('valid', 1));

// Sample your valid points for both PFSI and Soil Moisture values

var allSamples = ee.FeatureCollection([]); // to store samples 

// Loop through year-season pairs (e.g. Summer-2015)
validYearSeasonPairs.getInfo().forEach(function(pair) {
  var year = pair[0];
  var season = pair[1];
  var key = year + '_' + season;

  var pfsi = ee.Image(pfsi3x3Dict.get(key));
  var sm = ee.Image(sm3x3Dict.get(key));

  if (pfsi && sm) {
    var samplesForThisPair = filteredMasterSamples.map(function(f) {
      var geom = f.geometry();

      var pfsiVal = pfsi.reduceRegion({
        reducer: ee.Reducer.first(), 
        geometry: geom,
        scale: 30,
        maxPixels: 1e9
      }).get(key);

      var smVal = sm.reduceRegion({
        reducer: ee.Reducer.first(), 
        geometry: geom,
        scale: 30,
        maxPixels: 1e9
      }).get(key);

      return f.set({
        'year': year,
        'season': season,
        'PFSI': ee.Number(pfsiVal),
        'SM': smVal
      });
    });

    allSamples = allSamples.merge(samplesForThisPair);
  } 
}); 

// -----------------------------------------------------------------------------
// 5. EXPORT FINAL DATASET
// -----------------------------------------------------------------------------

Export.table.toDrive({
  collection: allSamples,
  description: 'Final_Dataset',
  fileFormat: 'CSV'
});

// -----------------------------------------------------------------------------
// 6. SAMPLE 2014 SOIL MOISTURE VALUES
// Vegetation response was expected to lag.
// -----------------------------------------------------------------------------

var sm2014Seasons = ['Summer','Autumn','Winter','Spring'];

var sm2014Samples = ee.FeatureCollection([]); // to store samples

// Loop over the four seasons
sm2014Seasons.forEach(function(season) {
  var key = '2014_' + season;

  // Grab the soil moisture image
  var sm3x3 = ee.Image(sm3x3Dict.get(key));

  // Sample
  var samplesForSeason = filteredMasterSamples.map(function(f) {
    var smVal = sm3x3.reduceRegion({
      reducer: ee.Reducer.first(),
      geometry: f.geometry(),
      scale: 30,
      maxPixels: 1e9
    }).get(key);      

    // Attach year, season, and the Soil Moisture value
    return f.set({
      'year':   2014,
      'season': season,
      'SM_2014': smVal
    });
  });

  sm2014Samples = sm2014Samples.merge(samplesForSeason);
});


// Export dataset 
Export.table.toDrive({
  collection: sm2014Samples,
  description: 'SoilMoisture_2014_SampledPoints',
  fileFormat: 'CSV'
});