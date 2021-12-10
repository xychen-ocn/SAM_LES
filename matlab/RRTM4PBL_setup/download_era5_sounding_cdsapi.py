#!/usr/bin/env python

import cdsapi


c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means',
    {
        'format': 'netcdf',
        'product_type': 'monthly_averaged_reanalysis',
        'variable': [
            'geopotential', 'ozone_mass_mixing_ratio', 'relative_humidity',
            'specific_humidity', 'temperature', 'u_component_of_wind',
            'v_component_of_wind', 'vertical_velocity',
        ],
        'pressure_level': [
            '50', '70', '100',
            '125', '150', '175',
            '200', '225', '250',
            '300', '350', '400',
            '450', '500', '550',
            '600', '650', '700',
            '750', '775', '800',
            '825', '850', '875',
            '900', '925', '950',
            '975', '1000',
        ],
        'year': '2004',
        'month': [
            '11', '12',
        ],
        'time': '00:00',
        'area': [
            20, -63, 15,
            57,
        ],
    },
    'era5_pressure_levels_monthly_means_rico_region_NovDec2004.nc')

