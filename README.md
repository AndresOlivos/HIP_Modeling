# Code for Habitat Intrinsic Potential Modeling of streams (for ArcGIS Pro 2.8)

### Description
This folder contains the code used for the modeling and analysis of stream habitats in drainage networks derived from Digital Elevation (DEM) and Precipitation models. We follow the hydrogeomorphic modeling procedures described in Clarke et al. (2008) for the later classification of reach habitats based on a Habitat Intrinsic Potential (HIP) model, originally proposed by Burnett et al. (2003, 2007) for mapping potential habitats for salmonids in Oregon. A detailed description of our objectives, methods, and justification of this approach will be provided in a peer-reviewed article Olivos et al. (2022): [DOI TBD], where we applied the HIP model to map potential habitats for non-native species in southern South America, including Chinook salmon, coho salmon, and North American beavers.

### Code
Most of these scripts operate through the user interface of ArcGIS Pro 2.X (ESRI, 2021), although individual scripts can also be used through any Python editing program, provided that ArcGIS Pro is installed. The use of these geoprocessing scripts require input data in metric units (DEM in meters above the level of the sea, and precipitation model in average annual cubic meters) and projected into a metric coordinate systems.

### Inputs
#### Layers
Digital Elevation Model
Digital Precipitation Model - Values representing average annual precipitation (cms)
Hydrological regions - Polygons representing the regional stratification for fitting dicharge formulas

#### Parameters
Minimum catchment area (for initial delineation)
Minimum annual discharge (for filtering initial channels)
Reach resolution (for line segmentation)
Parallel Processing Factor

### Outputs
#### Products
Streams (segmented polylines)
Points (point representation of reaches)

#### Sub-products
Filled DEM
Flow Direction raster
Flow Accumulation raster (unweighted)
Flow Accumulation raster (precipitation weighted)
Flow Accumulation raster (slope weighted)
Flow lines (unweighted, filtered)
Flow lines (precip. weighted - filtered)
Flow lines (slope weighted - filtered)
Line order (not strahler)

### Fields
