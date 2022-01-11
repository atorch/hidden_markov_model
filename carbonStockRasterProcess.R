##Process Carbon Stock tiff file and ensure 


##Carbon Stock File
carbonDir <- '/home/ted/Dropbox/amazon_hmm_shared/carbon_stock_data'
carbonStockRasterList <- lapply(list.files(carbonDir, pattern = 'W\\.tif$',full.names=TRUE),terra::rast)
carbonStockRaster1 <- do.call(terra::merge,carbonStockRasterList[1:3])
carbonStockRaster2 <- do.call(terra::merge,carbonStockRasterList[4:7])

carbonStockRaster <- terra::merge(carbonStockRaster1,
                                  carbonStockRaster2)

terra::writeRaster(carbonStockRaster, filename = file.path(carbonDir,'carbonStockRaster2000.tif'),tempdir = '')


carbon2017File <- '/home/ted/Dropbox/amazon_hmm_shared/CarbonStock_2017/pa_br_carbonMap_improved_ds_oskar_chalmers.tif'
carbonStockRaster2017Orig <- terra::rast(carbon2017File)
carbonStockRaster2017 <- terra::project(carbonStockRaster2017Orig,'+proj=longlat +datum=WGS84 +no_defs')
terra::writeRaster(carbonStockRaster2017, filename = file.path(carbonDir,'carbonStockRaster2017.tif'),tempdir = '')

