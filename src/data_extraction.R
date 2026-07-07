library(sf)

extract.tbe.data <- function(cfg) {
  tbe <- st_read(cfg$extraction$TBE$inpath)
  reg.adm <- st_read(cfg$extraction$SDA$inpath, layer = 'regio_s')
  
  region.of.interest <- subset(reg.adm, RES_NM_REG %in% cfg$extraction$SDA$regions) |>
    st_union()
  if (st_crs(region.of.interest) != st_crs(tbe)) {
    region.of.interest <- st_transform(region.of.interest, st_crs(tbe))
  }
  tbe.subset <- subset(tbe, st_intersects(tbe, region.of.interest, sparse = FALSE))
  
  st_write(tbe.subset, cfg$extraction$TBE$outpath, append = FALSE)
}


extract.ieqm.data <- function(cfg) {
  ieqm.list <- list.files(cfg$extraction$IEQM$inpath, full.names = TRUE)
  reg.adm <- st_read(cfg$extraction$SDA$inpath, layer = 'regio_s')
  
  region.of.interest <- subset(reg.adm, RES_NM_REG %in% cfg$extraction$SDA$regions) |>
    st_union()
  output.data <- NULL
  for(file in ieqm.list){
    ieqm <- st_read(file, layer = 'pee_ori')
    if (st_crs(region.of.interest) != st_crs(ieqm)) {
      region.of.interest <- st_transform(region.of.interest, st_crs(ieqm))
    }
    ieqm.subset <- subset(ieqm, st_intersects(region.of.interest, ieqm, sparse = FALSE))
    if(is.null(output.data)){
      output.data <- ieqm.subset
    } else {
      output.data <- rbind(output.data, ieqm.subset)
    }
  }
  output.data <- subset(output.data, !duplicated(output.data))
  st_write(output.data, cfg$extraction$IEQM$outpath, append = FALSE)
}
