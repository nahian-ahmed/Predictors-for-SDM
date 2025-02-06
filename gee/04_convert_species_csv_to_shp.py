import pandas as pd
import geopandas as gpd

if __name__ == "__main__":
    
    # Load in the species observations
    species_df = pd.read_csv("C:/Users/johnb/Downloads/OR2020_SpeciesRecords.csv")
    
    # Split the column by underscore, and take only latitude and longitude columns
    species_df[['Latitude','Longitude']] = (
        species_df['Unique_Checklist_ID']
        .str.split('_', expand=True)
        .iloc[:, [0, 1]]
        .astype(float)
        )
           
    # Convert to a GeoDataFrame
    species_gdf = gpd.GeoDataFrame(
        species_df,
        geometry = gpd.points_from_xy(species_df.Longitude, species_df.Latitude),
        crs = "EPSG:4326"
        )    
    
    # Subset just the ID and the geometry
    species_gdf = species_gdf[["Unique_Checklist_ID", "Year", "geometry"]]
    species_gdf = species_gdf.rename(columns={'Unique_Checklist_ID': 'Unique_ID'})
        
    # Save as a shapefile
    species_gdf.to_file("C:/Users/johnb/Downloads/OR2020_SpeciesRecords.shp", driver="ESRI Shapefile")        
    