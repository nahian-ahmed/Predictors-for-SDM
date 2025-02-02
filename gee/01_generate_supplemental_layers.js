/*************************************************************
*
* Title: 01_generate_supplemental_layers 
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

function main () {
  
    // Get the geometry of Oregon State
    var oregon = ee.FeatureCollection("TIGER/2016/States")
      .filter(ee.Filter.eq("NAME", "Oregon"))
      .geometry()
      .buffer(25000);
    
    // Create a binary mask
    var binary_mask = ee.Image(0).paint(oregon, 1);
    
    // Generate a covering grid
    var covering_grid = oregon.coveringGrid({
      proj: "EPSG:5070", 
      scale: 250000
    });
  
    // Export the covering grid
    var aoi_oregon = ee.FeatureCollection([
      ee.Feature(oregon, {"ID": "OREGON_25KM_BUFFER"})
    ]);
    Export.table.toAsset({
      collection: aoi_oregon, 
      description: "BirdSDM2025-Oregon-Oregon25kmBuffer",
      assetId: "users/JohnBKilbride/Bird_SDM_2025/oregon_25km_buffer" 
    });
    
    // Export the binary mask as an image
    Export.image.toAsset({
      image: binary_mask.byte(), 
      description: "BirdSDM2025-Oregon-Mask", 
      assetId: "users/JohnBKilbride/Bird_SDM_2025/oregon_mask", 
      region: oregon.bounds(), 
      scale: 30, 
      crs: "EPSG:5070",
      maxPixels: 1e13 
    });
    
    // Export the covering grid
    Export.table.toAsset({
      collection: covering_grid, 
      description: "BirdSDM2025-Oregon-250kmGrid",
      assetId: "users/JohnBKilbride/Bird_SDM_2025/oregon_250km_grid" 
    });
  
    return null;
  
  }
  
  main();
  