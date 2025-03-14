import copernicusmarine

copernicusmarine.subset(
  dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
  variables=["uo", "vo"],
  minimum_longitude=20,
  maximum_longitude=70,
  minimum_latitude=-40,
  maximum_latitude=10,
  start_datetime="2017-01-01T00:00:00",
  end_datetime="2021-06-30T00:00:00",
  minimum_depth=0,
  maximum_depth=1,
  output_filename = "CMEMS_Indian_currents_surface.nc",
  output_directory = "copernicus-data"
)