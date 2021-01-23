# 1. Create account for Copernicus Data Store
# 2. Install the CDS api key (save a .cdsapirc file with url and key in home directory)
# 3. Install the cds api on the computer (pip install cdsapi)

import cdsapi

c = cdsapi.Client()

import os

os.chdir('G:/Thesis_Data/Precipitation/')

fyear=1984  
lyear=2020   

   
for year in range(fyear,lyear):
    c.retrieve(
        'reanalysis-era5-land-monthly-means',
        {
            'format': 'netcdf',
            'variable': 'total_precipitation',
            'product_type': 'monthly_averaged_reanalysis',
            'area': '31.4/-7.5/31.0/-6.9',
            'year': str(year),
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                ],
            'time': '00:00',
            },
        str(year)+'.nc')
