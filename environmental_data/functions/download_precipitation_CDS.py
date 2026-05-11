# API request for data - Python:
import cdsapi

dataset = "reanalysis-era5-land-monthly-means"
request = {
  "product_type": ["monthly_averaged_reanalysis"],
  "variable": ["total_precipitation"],
  "year": [ "2004", "2005",
    "2006", "2007", "2008",
    "2009", "2010", "2011",
    "2012", "2013", "2014",
    "2015", "2016", "2017",
    "2018", "2019", "2020",
    "2021", "2022", "2023",
    "2024", "2025"
  ],
  "month": [
    "01", "02", "03",
    "04", "05", "06",
    "07", "08", "09",
    "10", "11", "12"
  ],
  "time": ["00:00"],
  "data_format": "GRIB",
  "area": [45.3831869151541, -11.8091501887745, 34.7325737729313, 4.24153981962996]
}
target = "environmental_data\data\CDS\precipitation.grib"


client = cdsapi.Client()
client.retrieve(dataset, request, target)

