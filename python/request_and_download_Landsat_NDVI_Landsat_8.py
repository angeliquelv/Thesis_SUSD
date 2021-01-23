from lsru import Usgs,Espa
from datetime import datetime

# Instantiate Espa and Usgs class 
espa=Espa()
usgs=Usgs()
usgs.login()

# Define extent
# Define query extent
bbox = (-7.45,31,-6.9,31.4)

collection = usgs.get_collection_name(8)

# Input for USGS api to find available scenes
meta_list_8 = usgs.search(collection=collection,
                         bbox=bbox,
                         begin=datetime(1984,1,1),
                         end=datetime(2020,1,1),
                         max_cloud_cover=100)

print(len(meta_list_8))

scene_list_8 = [x['displayId'] for x in meta_list_8]

### Place order to ESPA On-demand interface
order_8 = espa.order(scene_list_8,products=['sr_ndvi','pixel_qa'],format='gtiff')


# When order is complete, download the data
for order in espa.orders:
    if order.is_complete:
        order.download_all_complete('G:/Thesis_data/NDVI/Landsat8')




