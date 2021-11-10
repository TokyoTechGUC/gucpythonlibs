import gdal


def resample_tif(input, reference, output):
    gdal.UseExceptions()
    ds = None
    ds = gdal.Open(reference, gdal.GA_ReadOnly)
    pop = gdal.Open(reference, gdal.GA_ReadOnly)
    pop_prj = pop.GetProjection()
    pop_geotrans = pop.GetGeoTransform()
    wide = pop.RasterXSize
    high = pop.RasterYSize

    src = gdal.Open(input, gdal.GA_ReadOnly)
    src.GetRasterBand(1).SetNoDataValue(0.0)
    dst_filename = output
    DriverGTiff = gdal.GetDriverByName("GTiff")

    dst = DriverGTiff.Create(
        dst_filename,
        wide,
        high,
        1,
        gdal.GDT_Float32,
        options=["COMPRESS=LZW", "PREDICTOR=2", "BIGTIFF=YES"],
    )
    dst.SetGeoTransform(pop_geotrans)
    dst.SetProjection(pop_prj)
    dst.GetRasterBand(1).SetNoDataValue(0.0)
    prj = dst.GetProjection()
    gdal.ReprojectImage(src, dst, prj, pop_prj, gdal.GRA_Average)


def resample_sleuth_tif(input, reference, output):
    color_dict = {
        (249, 209, 110): 101.0,
        (20, 52, 214): 102.0,
        (255, 255, 255): 0.0,
        (60, 10, 0): 1.5,
        (60, 15, 0): 2.5,
        (60, 20, 0): 3.5,
        (60, 25, 0): 4.5,
        (60, 30, 0): 5.5,
        (60, 35, 0): 6.5,
        (60, 40, 0): 7.5,
        (60, 45, 0): 8.5,
        (60, 50, 0): 9.5,
        (60, 55, 0): 10.5,
        (60, 60, 0): 11.5,
        (60, 65, 0): 12.5,
        (60, 70, 0): 13.5,
        (60, 75, 0): 14.5,
        (60, 80, 0): 15.5,
        (60, 85, 0): 16.5,
        (60, 90, 0): 17.5,
        (60, 95, 0): 18.5,
        (60, 100, 0): 19.5,
        (60, 105, 0): 20.5,
        (60, 110, 0): 21.5,
        (60, 115, 0): 22.5,
        (60, 120, 0): 23.5,
        (60, 125, 0): 24.5,
        (60, 130, 0): 25.5,
        (60, 135, 0): 26.5,
        (60, 140, 0): 27.5,
        (60, 145, 0): 28.5,
        (60, 150, 0): 29.5,
        (60, 155, 0): 30.5,
        (60, 160, 0): 31.5,
        (60, 165, 0): 32.5,
        (60, 170, 0): 33.5,
        (60, 175, 0): 34.5,
        (60, 180, 0): 35.5,
        (60, 185, 0): 36.5,
        (60, 190, 0): 37.5,
        (60, 195, 0): 38.5,
        (60, 200, 0): 39.5,
        (60, 205, 0): 40.5,
        (60, 210, 0): 41.5,
        (60, 215, 0): 42.5,
        (60, 220, 0): 43.5,
        (60, 225, 0): 44.5,
        (60, 230, 0): 45.5,
        (60, 235, 0): 46.5,
        (60, 240, 0): 47.5,
        (60, 245, 0): 48.5,
        (60, 250, 0): 49.5,
        (0, 90, 0): 50.5,
        (0, 110, 0): 51.5,
        (0, 130, 0): 52.5,
        (0, 150, 0): 53.5,
        (0, 170, 0): 54.5,
        (0, 190, 0): 55.5,
        (0, 210, 0): 56.5,
        (0, 230, 0): 57.5,
        (0, 255, 0): 58.5,
        (0, 0, 90): 59.5,
        (0, 0, 110): 60.5,
        (0, 0, 130): 61.5,
        (0, 0, 150): 62.5,
        (0, 0, 170): 63.5,
        (0, 0, 190): 64.5,
        (0, 0, 210): 65.5,
        (0, 0, 230): 66.5,
        (0, 0, 255): 67.5,
        (0, 90, 60): 68.5,
        (0, 90, 90): 69.5,
        (0, 90, 120): 70.5,
        (0, 90, 150): 71.5,
        (0, 90, 180): 72.5,
        (0, 90, 210): 73.5,
        (0, 90, 255): 74.5,
        (0, 120, 60): 75.5,
        (0, 120, 90): 76.5,
        (0, 120, 120): 77.5,
        (0, 120, 150): 78.5,
        (0, 120, 180): 79.5,
        (0, 120, 210): 80.5,
        (0, 120, 255): 81.5,
        (90, 0, 30): 82.5,
        (90, 0, 60): 83.5,
        (90, 0, 90): 84.5,
        (90, 0, 120): 85.5,
        (90, 0, 150): 86.5,
        (90, 0, 180): 87.5,
        (90, 0, 210): 88.5,
        (90, 0, 255): 89.5,
        (150, 90, 0): 90.5,
        (150, 120, 0): 91.5,
        (150, 150, 0): 92.5,
        (150, 180, 0): 93.5,
        (150, 90, 255): 94.5,
        (255, 150, 0): 95.5,
        (255, 180, 0): 96.5,
        (255, 210, 0): 97.5,
        (255, 255, 0): 98.5,
        (255, 0, 51): 99.5,
    }
    import numpy as np
    from PIL import Image
    import numpy.lib.recfunctions as nlr

    gdal.UseExceptions()
    pop = gdal.Open(reference, gdal.GA_ReadOnly)
    pop_prj = pop.GetProjection()
    pop_geotrans = pop.GetGeoTransform()
    wide = pop.RasterXSize
    high = pop.RasterYSize

    im = Image.open(input)
    rgbimage = im.convert("RGB")
    data = np.asarray(rgbimage)
    data = nlr.unstructured_to_structured(data).astype("O")
    data = np.vectorize(color_dict.get)(data).astype("float")
    data = np.where(np.isnan(data), 0.0, data)

    dst_filename = output

    DriverGTiff = gdal.GetDriverByName("GTiff")
    dst = DriverGTiff.Create(
        dst_filename,
        wide,
        high,
        1,
        gdal.GDT_Float32,
        options=["COMPRESS=LZW", "PREDICTOR=2", "BIGTIFF=YES"],
    )
    dst.SetGeoTransform(pop_geotrans)
    dst.SetProjection(pop_prj)
    outband = dst.GetRasterBand(1)
    outband.WriteArray(data)


def copy_projection(input, reference, output):
    import numpy as np
    from PIL import Image
    import numpy.lib.recfunctions as nlr

    gdal.UseExceptions()
    pop = gdal.Open(reference, gdal.GA_ReadOnly)
    pop_prj = pop.GetProjection()
    pop_geotrans = pop.GetGeoTransform()
    wide = pop.RasterXSize
    high = pop.RasterYSize

    im = Image.open(input)
    # 	rgbimage=list(im.getdata())#convert('RGB')
    data = np.asarray(im)
    # data = nlr.unstructured_to_structured(data).astype('O')
    data = np.where(np.isnan(data), 0.0, data)
    print(data.shape)

    dst_filename = output

    DriverGTiff = gdal.GetDriverByName("GTiff")
    dst = DriverGTiff.Create(
        dst_filename,
        wide,
        high,
        1,
        gdal.GDT_Float32,
        options=["COMPRESS=LZW", "PREDICTOR=2", "BIGTIFF=YES"],
    )
    dst.SetGeoTransform(pop_geotrans)
    dst.SetProjection(pop_prj)
    outband = dst.GetRasterBand(1)
    outband.WriteArray(data)


def resample_sleuth_growth_tif(input, reference, output):
    color_dict = {
        (255, 0, 0): 1,  # Seed urban
        (0, 255, 0): 2,  # Diffusion
        (0, 0, 255): 6,  # Not used
        (255, 255, 0): 3,  # Breed
        (255, 255, 255): 4,  # Spread
        (0, 255, 255): 5,
        (0, 0, 0): 0,
    }  # Road-influenced

    import numpy as np
    from PIL import Image
    import numpy.lib.recfunctions as nlr

    gdal.UseExceptions()
    pop = gdal.Open(reference, gdal.GA_ReadOnly)
    pop_prj = pop.GetProjection()
    pop_geotrans = pop.GetGeoTransform()
    wide = pop.RasterXSize
    high = pop.RasterYSize

    im = Image.open(input)
    rgbimage = im.convert("RGB")
    data = np.asarray(rgbimage)
    data = nlr.unstructured_to_structured(data).astype("O")  # ;print(data);exit()
    data = np.vectorize(color_dict.get)(data).astype("int")
    data = np.where(np.isnan(data), 0.0, data)

    dst_filename = output

    DriverGTiff = gdal.GetDriverByName("GTiff")
    dst = DriverGTiff.Create(
        dst_filename,
        wide,
        high,
        1,
        gdal.GDT_Int16,
        options=["COMPRESS=LZW", "PREDICTOR=2", "BIGTIFF=YES"],
    )
    dst.SetGeoTransform(pop_geotrans)
    dst.SetProjection(pop_prj)
    outband = dst.GetRasterBand(1)
    outband.WriteArray(data)
