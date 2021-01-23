from lsru import Usgs,Espa
from datetime import datetime

# Instantiate Espa and Usgs class 
espa=Espa()
usgs=Usgs()
usgs.login()

# Define extent
# Define query extent
bbox = (-7.5,31,-7,31.4)

collection = usgs.get_collection_name(5)

# Input for USGS api to find available scenes
meta_list_5 = usgs.search(collection=collection,
                         bbox=bbox,
                         begin=datetime(1984,1,1),
                         end=datetime(2020,1,1),
                         max_cloud_cover=100)

print(len(meta_list_5))

scene_list_5 = [x['displayId'] for x in meta_list_5]

### Place order to ESPA On-demand interface
order_5 = espa.order(scene_list_5,products=['sr_ndvi','pixel_qa'],format='gtiff')


collection = usgs.get_collection_name(7)

# Input for USGS api to find available scenes
meta_list_7 = usgs.search(collection=collection,
                         bbox=bbox,
                         begin=datetime(1984,1,1),
                         end=datetime(2020,1,1),
                         max_cloud_cover=100)

print(len(meta_list_7))

scene_list_7 = [x['displayId'] for x in meta_list_7]

### Place order to ESPA On-demand interface
order_7 = espa.order(scene_list_7,products=['sr_ndvi','pixel_qa'],format='gtiff')

