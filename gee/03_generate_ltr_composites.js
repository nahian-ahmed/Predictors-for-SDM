/*************************************************************
*
* Title: 03_generate_ltr_composites
* Author: John Kilbride
* Date: 2025-02-01
* Description: 
*  This script generates LandTrendr interpolated composites.
* 
*************************************************************/

// Primary logic for the script
function main () {
  
  // Define the start year and the end year of the time series
  var start_year = 1990;
  var end_year = 2024;
  
  // Load in the AOI of Oregon w/ a 25km buffer
  var study_area = ee.FeatureCollection("users/JohnBKilbride/Bird_SDM_2025/oregon_25km_buffer").geometry();
  
  // // Create a temporary mask
  // var temp_mask = ee.Image(0).paint(geometry, 1);
  
  // Load in the BRDF-adjusted Medoids
  var medoids = ee.ImageCollection("users/JohnBKilbride/Bird_SDM_2025/oregon_brdf_medoids");
    // .map(function (image) {
    //   return image.updateMask(temp_mask);
    // });
  
  // Add the inverted NBR band for LandTrendr segmentation
  medoids = medoids.map(add_seg_band_nbr);
  
  // Run the LandTrendr interpolation
  var interpolated_medoids = interpolate_time_series(medoids, start_year, end_year);
  
  // Iterate over the years of interest
  for (var current_year = start_year; current_year <= end_year; current_year++) {
    
    // Get the export image
    var export_image = ee.Image(interpolated_medoids.filterDate(current_year+"-01-01", current_year+"-12-31").first());
      
    // Create the metadata dates
    var start_date = ee.Date.fromYMD(current_year, 1, 1).millis();
    var end_date = ee.Date.fromYMD(current_year, 12, 31).millis();
    
    // Assign the metadata 
    export_image = export_image.set({
      "system:time_start": start_date,
      "system:time_end": end_date
    });
    
    // Export the binary mask as an image
    Export.image.toAsset({
      image: export_image.toInt16(), 
      description: "BirdSDM2025-Oregon-MedoidBRDF-" + current_year, 
      assetId: "users/JohnBKilbride/Bird_SDM_2025/oregon_interpolated_medoids/interpolated_" + current_year, 
      region: study_area.bounds(), 
      scale: 30, 
      crs: "EPSG:5070",
      maxPixels: 1e13 
    });

  }
  
  // var viz =  {min:0, max:[1500,5500,1500], bands:["B7","B4","B3"]};
  // Map.addLayer(medoids, viz, "Medoids");
  // Map.addLayer(interpolated_medoids, viz, "Interpolated Medoids");
  // Map.addLayer(export_image, viz, "Medoids");


  return null;
  
}

// Add an inverted NBR band for segmentation
function add_seg_band_nbr(image) {

  // Compute the inverted NBR band
  var nbr = image.normalizedDifference(['B4', 'B7'])
    .rename('NBR_SEG')
    .multiply(-1);
  
  return nbr.addBands(image).set({
    'system:time_start': image.date().millis()
    }).float();

}

// Run the LandTrendr algorithm and performs the interpolation logic
function interpolate_time_series (image_collection, start_year, end_year) {

  // Run the LandTrendr algorithm
  var ltr_output = ee.Algorithms.TemporalSegmentation.LandTrendr({
    timeSeries: image_collection,
    maxSegments: 10,
    spikeThreshold: 0.9,
    vertexCountOvershoot: 3,
    preventOneYearRecovery: true,
    recoveryThreshold: 0.25,
    pvalThreshold: 0.05,
    bestModelProportion: 0.75,
    minObservationsNeeded: 6
  });

  // Flatten the LandTrendr outputs
  var flattened_images = flatten_ltr_fits(ltr_output, start_year, end_year);

  return flattened_images;
  
}

// Flatten the LandTrendr outputs
function flatten_ltr_fits (ltr_output, start_year, end_year) {

  // Cast the years as ee.Numbers to get the index values
  var ee_start_year = ee.Number(start_year);
  var ee_end_year = ee.Number(end_year);

  // Compute the start and end indices to slice the LandTrendr outputs
  var start_index = 0;
  var end_index = ee_end_year.subtract(ee_start_year);

  // Select the fitted bands of the LandTrendr output
  var fitted_values = ltr_output
      .select(['B1_fit', 'B2_fit', 'B3_fit', 'B4_fit', 'B5_fit', 'B7_fit'])
      .toArray(1);

  // Define an ee.List of the years
  var years = ee.List.sequence(start_index, end_index);

  // Function to flatten a single year from the LandTrendr fitted values
  var flatten_single_year = function(index) {
    
    // Cast the input year as an ee.Number and define the current date
    index = ee.Number(index);
    var date = ee.Date.fromYMD(index.add(ee_start_year), 7, 1).millis();

    // Define the layer used to unmask the fitted values
    var unmask_layer = ee.Image([0, 0, 0, 0, 0, 0])
        .toArray()
        .toArray(1)
        .arrayTranspose();

    // Slice out the values of the given year
    var year_fit = fitted_values
        .unmask(unmask_layer)
        .arraySlice(0, index.toInt16(), index.add(1).toInt16())
        .arrayProject([1]);

    // Flatten the fitted spectral values for a single year
    return year_fit
        .arrayFlatten([['B1', 'B2', 'B3', 'B4', 'B5', 'B7']])
        .unmask(ee.Image([0, 0, 0, 0, 0, 0]))
        .set('system:time_start', date)
        .toFloat();
};

  // Map the flattening function over the list of years to get a 
  // time-series of interpolated images
  var outputs = ee.ImageCollection.fromImages(years.map(flatten_single_year));

  return outputs;

}


main();

