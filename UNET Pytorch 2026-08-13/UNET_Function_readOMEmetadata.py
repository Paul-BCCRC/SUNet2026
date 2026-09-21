#Writen by Paul Gallagher (BC Cancer Research Centre)

def readOMEmetadata(ome, TIFF_page_number):
    print(r"Reading OME Metadata...")
    try:
        Objective = ome.instruments[0].objectives[0].calibrated_magnification
    except:
        Objective = None
    try:
        Resolution = ome.images[0].pixels.physical_size_x
    except:
        Resolution = None

    try:
        ImageType = ome.datasets[0].name
    except:
        ImageType = None

    try:
        ScaleF = ome.images[0].pixels.channels[TIFF_page_number-1].nd_filter
    except:
        ScaleF = None

    try:
        #print(ome.images[0].pixels)
        #print(ome.images[0].pixels.channels)
        N = len(ome.images[0].pixels.planes)
        ImageLabels = []
        for n in range(0, N):
            ImageLabels.append(ome.images[0].pixels.channels[n].name)
    except:
        n=1#just do something

    return [Objective, Resolution, ImageType, ScaleF, ImageLabels]