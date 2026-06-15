library(imageFeatureTCGA)

# 2. Fetch your URLs
url_h5ad    <- getCatalog("hovernet", format = "h5ad")    |> getFileURLs()
url_geojson <- getCatalog("hovernet", format = "geojson") |> getFileURLs()
url_thumb   <- getCatalog("hovernet", format = "thumb")   |> getFileURLs()
url_json    <- getCatalog("hovernet", format = "json")    |> getFileURLs()

saveRDS(url_h5ad, file = '~/ppp_TCGA/data/catalogs/h5ad_TCGAcatalog.rds') 
saveRDS(url_geojson, file = '~/ppp_TCGA/data/catalogs/geojson_TCGAcatalog.rds') 
saveRDS(url_thumb, file = '~/ppp_TCGA/data/catalogs/thumb_TCGAcatalog.rds') 
saveRDS(url_json, file = '~/ppp_TCGA/data/catalogs/json_TCGAcatalog.rds') 
