/*************************************************************
*
* Title: 02_generate_brdf_composites 
* Author: John Kilbride
* Date: 2025-01-29
* Description: 
*  This script generates a number of data products used in the creation of the
*  BRDF-adjusted Landsat satellite images including:
*    - A mask of the areas to process
*    - A set of 150 km that will be processed
*    - 
* 
*************************************************************/

// Load in the BRDF logic
var BRDF = require("users/JohnBKilbride/Bird_SDM_2025/:module/brdf_logic");

// Primary compositing logic
function main () {
  
  // Load in the AOI of Oregon w/ a 25km buffer
  var study_area = ee.FeatureCollection("users/JohnBKilbride/Bird_SDM_2025/oregon_25km_buffer").geometry();
  
  // Load in the binary mask of Oregon
  var binary_mask = ee.Image("users/JohnBKilbride/Bird_SDM_2025/oregon_mask");
  
  // Define the start year and the end year
  var start_year = 1990;
  var end_year = 2024;
  
  // Define the start day and the end day of the compositing period
  var start_day = 183;
  var end_day = 244;
  
  // Iterate over the years of interest
  for (var current_year = start_year; current_year <= end_year; current_year++) {
    
    // Load in the data for the medoid composite
    var annual_ls_collection = load_landsat_imagery(current_year, study_area, start_day, end_day);
    
    // Compute the medoid
    var medoid_composite = ee.Image(make_medoid_composite(annual_ls_collection))
      .updateMask(binary_mask)
      .toInt16();
      
    // Create the metadata dates
    var start_date = ee.Date.fromYMD(current_year, 1, 1).advance(start_day, 'day').millis();
    var end_date = ee.Date.fromYMD(current_year, 1, 1).advance(end_day, 'day').millis();
    
    // Assign the metadata 
    medoid_composite = medoid_composite.set({
      "system:time_start": start_date,
      "system:time_end": end_date
    });
    
    // Export the binary mask as an image
    Export.image.toAsset({
      image: medoid_composite, 
      description: "BirdSDM2025-Oregon-MedoidBRDF-" + current_year, 
      assetId: "users/JohnBKilbride/Bird_SDM_2025/oregon_brdf_medoids/medoid_" + current_year, 
      region: study_area.bounds(), 
      scale: 30, 
      crs: "EPSG:5070",
      maxPixels: 1e13 
    });
  
  }
  
  Map.addLayer(medoid_composite, {min:0, max:[1500,5500,1500], bands:['B7','B4','B3']}, "Medoid " + current_year);

  
  return null;
  
}

// Load the imagery for a given period of time
function load_landsat_imagery (composite_year, study_area_geo, start_day, end_day) {

  // Construct the ee.Dates to be used when filtering the time-series
  var start_date = ee.Date.fromYMD(composite_year, 1, 1);
  var end_date = ee.Date.fromYMD(composite_year, 12, 31);
  
  // Load the Landsat Imagery
  var landsat_5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
    .filterDate(start_date, end_date)
    .filter(ee.Filter.calendarRange(start_day, end_day))
    .filterBounds(study_area_geo)
    .map(apply_scale_factors)
    .select(
        ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7", "QA_PIXEL"],
        ["B1", "B2", "B3", "B4", "B5", "B7", "QA_PIXEL"]
    )
    .map(BRDF.apply_brdf)
    .map(apply_pixel_qa);
  
  var landsat_7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
    .filterDate(start_date, end_date)
    .filter(ee.Filter.calendarRange(start_day, end_day))
    .filterBounds(study_area_geo)
    .map(apply_scale_factors)
    .select(
        ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7", "QA_PIXEL"],
        ["B1", "B2", "B3", "B4", "B5", "B7", "QA_PIXEL"]
    )
    .map(BRDF.apply_brdf)
    .map(apply_pixel_qa);
      
  var landsat_8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
    .filterDate(start_date, end_date)
    .filter(ee.Filter.calendarRange(start_day, end_day))
    .filterBounds(study_area_geo)
    .map(apply_scale_factors)
    .select(
        ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "QA_PIXEL"],
        ["B1", "B2", "B3", "B4", "B5", "B7", "QA_PIXEL"]
    )
    .map(BRDF.apply_brdf)
    .map(apply_pixel_qa);
  
  var landsat_9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
    .filterDate(start_date, end_date)
    .filter(ee.Filter.calendarRange(start_day, end_day))
    .filterBounds(study_area_geo)
    .map(apply_scale_factors)
    .select(
        ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "QA_PIXEL"],
        ["B1", "B2", "B3", "B4", "B5", "B7", "QA_PIXEL"]
    )
    .map(BRDF.apply_brdf)
    .map(apply_pixel_qa);

  return ee.ImageCollection(landsat_9.merge(landsat_8).merge(landsat_7).merge(landsat_5));

}

// Apply the Landsat collection 2 scaling factors then scale by 10,000
function apply_scale_factors (image) {
  
  var spectral_bands = image.select('SR_B.*')
    .multiply(0.0000275)
    .add(-0.2)
    .multiply(10000);
  
  return image.addBands(spectral_bands, null, true).toInt16();

}

// Apply the Pixel QA band
function apply_pixel_qa(image) {

  // Extract the QA mask and apply the bitwise operators
  var qa = image.select("QA_PIXEL");
  var bit_mask = qa.bitwiseAnd(1 << 2).eq(0)  // Cirrus  
    .and(qa.bitwiseAnd(1 << 3).eq(0))         // Clouds
    .and(qa.bitwiseAnd(1 << 4).eq(0))         // Cloud shadow
    .and(qa.bitwiseAnd(1 << 5).eq(0));        // Snow

  // Drop the Pixel QA band
  image = image.select('B.*');

  return image.updateMask(bit_mask);

}

// Define the Medoid compositing function
function make_medoid_composite (input_collection) {
  
  // Add a dummy image to prevent issues with an empty collection
  var dummy_collection = ee.ImageCollection([
    ee.Image([0, 0, 0, 0, 0, 0]).mask(ee.Image(0))
  ]).select([0, 1, 2, 3, 4, 5], ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']);
  var images_to_composite = input_collection.merge(dummy_collection);

  // Compute the median across images in the collection per band
  var median_image = images_to_composite.median();

  // Calculate the distance of each observation from the median
  function inner_median_diff(img) {
    var difference = img.subtract(median_image).pow(2);
    return difference.reduce(ee.Reducer.sum()).addBands(img);
  }

  var median_difference = images_to_composite.map(inner_median_diff);

  // Select the image with the smallest difference from the median
  return ee.ImageCollection(median_difference)
    .reduce(ee.Reducer.min(7))
    .select([1, 2, 3, 4, 5, 6], ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']);
}





main();
