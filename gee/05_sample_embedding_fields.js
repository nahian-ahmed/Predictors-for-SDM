/*************************************************************
*
* Title: 03_generate_ltr_composites
* Author: John Kilbride
* Date: 2025-02-02
* Description: 
*  This script intersects the sample points w/ the embedding fields.
* 
*************************************************************/

function main () {
  
    // Define the study area -- Oregon
    var study_area = ee.FeatureCollection('TIGER/2018/States')
      .filter(ee.Filter.eq('NAME', 'Oregon'))
      .geometry();
    
    // Load in the embedding layers
    var embedding_collection = ee.ImageCollection('projects/mldp-partners/assets/preview/efm_v2_preview')
      .filterBounds(study_area);
  
    // Load in the species feature collection
    var species_data = ee.FeatureCollection('users/JohnBKilbride/Bird_SDM_2025/OR2020_SpeciesRecords');
    
    // Get the years to process
    var years = species_data.aggregate_array('Year')
      .distinct()
      .sort()
      .getInfo();
    
    // Iterate over the years to sample
    var output_featcol = ee.FeatureCollection([]);
    for (var i = 0; i < years.length; i++) {
      
      var species_year = years[i];
  
      // Set years prior to 2017 to 2017
      var embedding_year = species_year < 2017 ? 2017 : species_year;
      
      // Get the current embedding year
      var start_date = ee.Date.fromYMD(embedding_year, 1, 1);
      var end_date = ee.Date.fromYMD(embedding_year, 12, 31);
      var embedding_image = ee.Image(embedding_collection.filterDate(start_date, end_date).mosaic())
        .unmask(-9999);
      
      // // Get the species observations for the current year
      var species_data_subset = species_data.filter(ee.Filter.eq('Year', species_year));
  
      // Run the sampling
      var samples = embedding_image.sampleRegions({
          collection: species_data_subset,
          scale: 10,
          projection: 'EPSG:4326',
          geometries: false
      });
  
      output_featcol = output_featcol.merge(samples);
      
    }
    
    // Exporting to Drive is not available in JavaScript API, so print instead
    Export.table.toDrive({
      collection: output_featcol, 
      description: "Sampled-Embeddings", 
      fileNamePrefix: "OR2020_SpeciesRecords_EmbeddingFields", 
      fileFormat: "CSV"
    });
  
    return null;
  
  }
  
  main();
  