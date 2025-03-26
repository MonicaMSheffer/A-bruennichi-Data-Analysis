###Function for plotting maps of bioclimatic variables

myplot <- function(bio_variable_to_plot, df){
  
  mycontinents2=dplyr::filter(world,continent %in% c("Asia","Europe"))
  
  df_to_plot <- spider_final[,c("population_name", bio_variable_to_plot, "geom")]  
  df_spat <- st_as_sf(df_to_plot)
  
  p<- tm_shape(mycontinents2,bbox=tmaptools::bb(matrix(c(-5,40,45,60),2,2))) +
    tm_polygons() +
    tm_shape(df_spat) +
    tm_symbols(col = bio_variable_to_plot, fill = bio_variable_to_plot, size=.5,palette="RdYlGn") 
  
  return(print(p))
  
}
