# ==================================
# Author: Jack Kilbride, Tristan Goodbody & Elias Ayrey - Renoster Systems Inc.
# Date: 2024-01-07
# Description: landsat processing
# ==================================
import math

import ee

# This is a port of the JavaScript GEE implementation of the BRDF correction proceedure
# BRDF Definitions -- from https://code.earthengine.google.com/3a6761dea6f1bf54b03de1b84dc375c6


def apply_brdf_correction_to_image(image):
    """
    Applies a Bidirectional Reflectance Distribution Function (BRDF) correction to an
    input image. This function selects specific bands, computes view and solar angles,
    and adjusts the bands based on the BRDF coefficients.

    Args:
        image: An Earth Engine image object that contains the bands to be corrected.

    Returns:
        An Earth Engine image object with BRDF-corrected bands.

    Examples:
        # Assume `eeImage` is an Earth Engine image with required bands
        corrected_image = apply_brdf_correction_to_image(eeImage)
    """
    # Select specific bands and rename them for BRDF correction
    image = image.select(
        ["B1", "B2", "B3", "B4", "B5", "B7", "pixel_qa"],
        ["blue", "green", "red", "nir", "swir1", "swir2", "pixel_qa"],
    )

    # BRDF coefficients for different bands
    coefficientsByBand = {
        "blue": {"fiso": 0.0774, "fgeo": 0.0079, "fvol": 0.0372},
        "green": {"fiso": 0.1306, "fgeo": 0.0178, "fvol": 0.0580},
        "red": {"fiso": 0.1690, "fgeo": 0.0227, "fvol": 0.0574},
        "nir": {"fiso": 0.3093, "fgeo": 0.0330, "fvol": 0.1535},
        "swir1": {"fiso": 0.3430, "fgeo": 0.0453, "fvol": 0.1154},
        "swir2": {"fiso": 0.2658, "fgeo": 0.0387, "fvol": 0.0639},
    }

    # Calculate image geometry and solar position parameters
    corners = findCorners(image)
    viewAngles(image, corners)
    solarPosition(image)
    sunZenOut(image)

    # Compute Ross-Thick and Li-Thin kernels for BRDF correction
    set_p(image, "relativeSunViewAz", "i.sunAz - i.viewAz")
    rossThick(image, "kvol", "i.sunZen", "i.viewZen", "i.relativeSunViewAz")
    rossThick(image, "kvol0", "i.sunZenOut", 0, 0)
    liThin(image, "kgeo", "i.sunZen", "i.viewZen", "i.relativeSunViewAz")
    liThin(image, "kgeo0", "i.sunZenOut", 0, 0)

    # Adjust band values based on BRDF coefficients
    adjustBands(image, coefficientsByBand)

    # Revert the band names to their original names
    return image.select(
        ["blue", "green", "red", "nir", "swir1", "swir2", "pixel_qa"],
        ["B1", "B2", "B3", "B4", "B5", "B7", "pixel_qa"],
    )


def to_image(value):
    """
    Converts a given value into a constant Earth Engine image.

    Args:
        value: The value to be converted into an Earth Engine image.

    Returns:
        An Earth Engine image object with the given value as a constant value across all pixels.

    Examples:
        # Convert a numeric value to an Earth Engine image
        image = to_image(42)
    """
    constant_image = ee.Image.constant(value)
    return constant_image


def set_p(image, name, to_add, args=None):
    """
    Adds a new band to an Earth Engine image, optionally using provided arguments
    for the band calculation.

    Args:
        image: The Earth Engine image to which the band is to be added.
        name: The name of the new band.
        to_add: The value or expression to be used for the new band. If `args` is not None,
                `to_add` is assumed to be an expression using `args`.
        args: Optional. Additional arguments used in the expression for the new band.

    Returns:
        An Earth Engine image with the new band added.

    Examples:
        # Add a constant band to an image
        new_image = set_p(eeImage, 'new_band', 10)

        # Add a band based on an expression
        new_image = set_p(eeImage, 'diff_band', 'b1 - b2', args=['b1', 'b2'])
    """
    if args is not None:
        # Convert the expression to an image using the provided arguments
        to_add = toImage(image, to_add, args)
    else:
        # Convert the value to an image
        to_add = to_image(to_add)
    # Add the new band to the image
    image = image.addBands(to_add.rename(name), None, True)
    return image


def set_pIf(image, name, condition, trueValue, falseValue=None):
    """
    Conditionally adds a new band to an Earth Engine image based on a given condition.
    The new band is set to `trueValue` where the condition is true, and `falseValue`
    otherwise.

    Args:
        image: The Earth Engine image to which the conditional band is to be added.
        name: The name of the new band.
        condition: The condition used to determine the band's value. If true, `trueValue`
                   is used; otherwise, `falseValue` is used.
        trueValue: The value to set for pixels where the condition is true.
        falseValue: Optional. The value to set for pixels where the condition is false.
                    If not provided, defaults to None.

    Examples:
        # Add a conditional band to an image
        new_image = set_pIf(eeImage, 'new_band', 'b1 > 10', 1, 0)
    """
    # Convert condition, trueValue, and falseValue to images
    condition = to_image(condition)
    trueMasked = to_image(trueValue).mask(condition)
    falseMasked = to_image(falseValue).mask(invertMask(condition))

    # Combine true and false values based on the condition
    value = trueMasked.unmask(falseMasked)

    # Add the conditional band to the image
    set_p(image, name, value)


def invertMask(mask):
    """
    Inverts a mask, changing 1s to 0s and vice versa.

    Args:
        mask: An Earth Engine image representing the mask to be inverted.

    Returns:
        An Earth Engine image representing the inverted mask.

    Examples:
        # Invert a mask
        inverted_mask = invertMask(maskImage)
    """
    # Invert the mask values
    return mask.multiply(-1).add(1)


def toImage(image, band, args):
    """
    Converts a band or expression to an Earth Engine image.

    Args:
        image: The Earth Engine image from which the band or expression is derived.
        band: A string representing the band or expression to be converted.
        args: A list of arguments used in the expression.

    Returns:
        An Earth Engine image derived from the specified band or expression.

    Examples:
        # Select a band from an image
        band_image = toImage(eeImage, 'B1', [])

        # Use an expression to create an image
        expr_image = toImage(eeImage, 'b1 * 2', ['b1'])
    """
    if isinstance(band, str):
        # Evaluate expression or select band
        if "." in band or " " in band or "{" in band:
            band = image.expression(
                format_string({"pi": math.pi}, band, args), {"i": image}
            )
        else:
            band = image.select(band)
    return ee.Image(band)


def format_string(constants, s, args=None):
    """
    Formats a string using given constants and additional arguments.

    Args:
        constants: A dictionary of constant values used in the string formatting.
        s: The string to be formatted.
        args: Optional. A dictionary of additional arguments to be used in formatting.
              If None, an empty dictionary is used.

    Returns:
        A formatted string with constants and arguments replaced. If an undeclared
        argument is found, None is returned and an error message is printed.

    Examples:
        # Example with constants and arguments
        formatted = format_string({'pi': 3.14}, 'Area of circle with radius {r} is {pi}*{r}**2', {'r': 2})
    """
    if args is None:
        args = {}
    # Merge constants and arguments
    all_args = {**constants, **args}

    result = s
    while "{" in result:
        result = result.format_map(all_args)
        unresolved_args = [arg[1:-1] for arg in result.split("{")[1:]]
        for arg in unresolved_args:
            if arg not in all_args:
                print(f"Undeclared argument: {arg} s: {s} args: {args}")
                return None
    return result


def merge(o1, o2):
    """
    Merges two dictionaries into a single dictionary.

    Args:
        o1: The first dictionary to merge.
        o2: The second dictionary to merge.

    Returns:
        A merged dictionary containing keys and values from both input dictionaries.

    Examples:
        # Merge two dictionaries
        merged_dict = merge({'a': 1, 'b': 2}, {'c': 3})
    """

    def addAll(target, toAdd):
        # Add keys and values from toAdd to target
        for key in toAdd:
            target[key] = toAdd[key]

    result = {}
    addAll(result, o1)
    addAll(result, o2)
    return result


def findCorners(image):
    """
    Finds the corner coordinates of the bounding box of a given Earth Engine image.

    Args:
        image: The Earth Engine image whose corner coordinates are to be found.

    Returns:
        A dictionary with keys 'upperLeft', 'upperRight', 'lowerRight', 'lowerLeft',
        representing the corner coordinates of the image's bounding box.

    Examples:
        # Find corners of an Earth Engine image
        corners = findCorners(eeImage)
    """
    # Get image footprint and bounds
    footprint = ee.Geometry(image.get("system:footprint"))
    bounds = ee.List(footprint.bounds().coordinates().get(0))
    coords = footprint.coordinates()

    # Extract x and y coordinates
    xs = coords.map(lambda item: x(item))
    ys = coords.map(lambda item: y(item))

    def findCorner(targetValue, values):
        # Find the corner closest to the target value
        diff = values.map(
            lambda value: ee.Number(value).subtract(targetValue).abs()
        )
        minValue = diff.reduce(ee.Reducer.min())
        idx = diff.indexOf(minValue)
        return coords.get(idx)

    # Identify each corner based on its position
    lowerLeft = findCorner(x(bounds.get(0)), xs)
    lowerRight = findCorner(y(bounds.get(1)), ys)
    upperRight = findCorner(x(bounds.get(2)), xs)
    upperLeft = findCorner(y(bounds.get(3)), ys)
    return {
        "upperLeft": upperLeft,
        "upperRight": upperRight,
        "lowerRight": lowerRight,
        "lowerLeft": lowerLeft,
    }


def x(point):
    """
    Extracts the x-coordinate from a point.

    Args:
        point: The point from which the x-coordinate is extracted.

    Returns:
        The x-coordinate of the given point.

    Examples:
        # Get x-coordinate from a point
        x_coord = x([3, 4])
    """
    return ee.Number(ee.List(point).get(0))


def y(point):
    """
    Extracts the y-coordinate from a point.

    Args:
        point: The point from which the y-coordinate is extracted.

    Returns:
        The y-coordinate of the given point.

    Examples:
        # Get y-coordinate from a point
        y_coord = y([3, 4])
    """
    return ee.Number(ee.List(point).get(1))


def pointBetween(pointA, pointB):
    """
    Calculates the centroid of the line connecting two points.

    Args:
        pointA: The first point as a list or tuple of coordinates [x, y].
        pointB: The second point as a list or tuple of coordinates [x, y].

    Returns:
        The coordinates of the centroid as a list [x, y].

    Examples:
        # Find the midpoint between two points
        midpoint = pointBetween([1, 2], [3, 4])
    """
    # Create a line and find its centroid
    return ee.Geometry.LineString([pointA, pointB]).centroid().coordinates()


def slopeBetween(pointA, pointB):
    """
    Calculates the slope of the line between two points.

    Args:
        pointA: The first point as a list or tuple of coordinates [x, y].
        pointB: The second point as a list or tuple of coordinates [x, y].

    Returns:
        The slope of the line.

    Examples:
        # Calculate the slope between two points
        slope = slopeBetween([1, 2], [3, 4])
    """
    # Calculate the slope between pointA and pointB
    return (y(pointA).subtract(y(pointB))).divide(
        (x(pointA)).subtract(x(pointB))
    )


def toLine(pointA, pointB):
    """
    Creates a line string from two points.

    Args:
        pointA: The first point as a list or tuple of coordinates [x, y].
        pointB: The second point as a list or tuple of coordinates [x, y].

    Returns:
        An Earth Engine LineString object representing the line.

    Examples:
        # Create a LineString from two points
        line = toLine([1, 2], [3, 4])
    """
    # Create a LineString geometry from pointA and pointB
    return ee.Geometry.LineString([pointA, pointB])


def liThin(image, bandName, sunZen, viewZen, relativeSunViewAz):
    """
    Applies the Li-Thin kernel to an Earth Engine image for a specified band.

    Args:
        image: The Earth Engine image to be processed.
        bandName: The name of the band to which the kernel is applied.
        sunZen: The name of the band representing the solar zenith angle.
        viewZen: The name of the band representing the view zenith angle.
        relativeSunViewAz: The name of the band representing the relative azimuth
                           angle between the sun and view directions.

    Examples:
        # Apply Li-Thin kernel to an Earth Engine image
        liThin(eeImage, 'kernel_band', 'sunZen', 'viewZen', 'relAzimuth')
    """
    # Define arguments for the Li-Thin kernel
    args = {
        "sunZen": sunZen,
        "viewZen": viewZen,
        "relativeSunViewAz": relativeSunViewAz,
        "h/b": 2,
    }

    # Calculate auxiliary angles and distances for the Li-Thin kernel
    anglePrime(image, "sunZenPrime", sunZen)
    anglePrime(image, "viewZenPrime", viewZen)
    cosPhaseAngle(
        image,
        "cosPhaseAngle",
        "i.sunZenPrime",
        "i.viewZenPrime",
        relativeSunViewAz,
    )
    set_p(
        image,
        "distance",
        "sqrt(pow(tan(i.sunZenPrime), 2) + pow(tan(i.viewZenPrime), 2)"
        + "- 2 * tan(i.sunZenPrime) * tan(i.viewZenPrime) * cos({relativeSunViewAz}))",
        args,
    )
    set_p(image, "temp", "1/cos(i.sunZenPrime) + 1/cos(i.viewZenPrime)")
    set_p(
        image,
        "cosT",
        toImage(
            image,
            "{h/b} * sqrt(pow(i.distance, 2) + pow(tan(i.sunZenPrime) * tan(i.viewZenPrime) * sin({relativeSunViewAz}), 2))"
            + "/ i.temp",
            args,
        ).clamp(-1, 1),
    )
    set_p(image, "t", "acos(i.cosT)")
    set_p(image, "overlap", "(1/{pi}) * (i.t - sin(i.t) * i.cosT) * (i.temp)")
    set_pIf(image, "overlap", "i.overlap > 0", 0)
    set_p(
        image,
        bandName,
        "i.overlap - i.temp"
        + "+ (1/2) * (1 + i.cosPhaseAngle) * (1/cos(i.sunZenPrime)) * (1/cos(i.viewZenPrime))",
    )


def anglePrime(image, name, angle):
    """
    Computes the angle prime for a given angle in an Earth Engine image.

    Args:
        image: The Earth Engine image to which the computation is applied.
        name: The name of the band to store the computed angle prime.
        angle: The name of the band representing the angle for which the angle prime is computed.

    Examples:
        # Compute angle prime for a band in an Earth Engine image
        anglePrime(eeImage, 'anglePrimeBand', 'angleBand')
    """
    # Define arguments for the computation
    args = {"b/r": 1, "angle": angle}

    # Compute tangent of the angle prime and add it as a band
    set_p(image, "tanAnglePrime", "{b/r} * tan({angle})", args)

    # Set negative values of tanAnglePrime to 0
    set_pIf(image, "tanAnglePrime", "i.tanAnglePrime < 0", 0)

    # Compute and add angle prime as a band
    set_p(image, name, "atan(i.tanAnglePrime)")


def cosPhaseAngle(image, name, sunZen, viewZen, relativeSunViewAz):
    """
    Computes the cosine of the phase angle for specified bands in an Earth Engine image.

    Args:
        image: The Earth Engine image to which the computation is applied.
        name: The name of the band to store the computed cosine phase angle.
        sunZen: The name of the band representing the solar zenith angle.
        viewZen: The name of the band representing the view zenith angle.
        relativeSunViewAz: The name of the band representing the relative azimuth angle.

    Examples:
        # Compute cosine phase angle for specified bands in an Earth Engine image
        cosPhaseAngle(eeImage, 'cosPhaseAngleBand', 'sunZen', 'viewZen', 'relAzimuth')
    """
    # Define arguments for the computation
    args = {
        "sunZen": sunZen,
        "viewZen": viewZen,
        "relativeSunViewAz": relativeSunViewAz,
    }

    # Compute and clamp the cosine of the phase angle between -1 and 1
    set_p(
        image,
        name,
        toImage(
            image,
            "cos({sunZen}) * cos({viewZen})"
            + "+ sin({sunZen}) * sin({viewZen}) * cos({relativeSunViewAz})",
            args,
        ).clamp(-1, 1),
    )


def adjustBands(image, coefficientsByBand):
    """
    Adjusts the bands of an Earth Engine image using BRDF coefficients.

    Args:
        image: The Earth Engine image to be adjusted.
        coefficientsByBand: A dictionary containing BRDF coefficients for each band.
    """
    # Apply BRDF coefficients to each band
    for bandName in coefficientsByBand:
        applyCFactor(image, bandName, coefficientsByBand[bandName])


def applyCFactor(image, bandName, coefficients):
    """
    Applies a correction factor to a band of an Earth Engine image based on BRDF coefficients.

    Args:
        image: The Earth Engine image to be processed.
        bandName: The name of the band to which the correction factor is applied.
        coefficients: A dictionary of BRDF coefficients for the specified band.

    Examples:
        # Apply correction factor to a band
        coefficients = {'fiso': 0.0774, 'fgeo': 0.0079, 'fvol': 0.0372}
        applyCFactor(eeImage, 'blue', coefficients)
    """
    # Compute BRDF and BRDF0 values
    brdf(image, "brdf", "kvol", "kgeo", coefficients)
    brdf(image, "brdf0", "kvol0", "kgeo0", coefficients)

    # Calculate and apply the correction factor
    set_p(image, "cFactor", "i.brdf0 / i.brdf", coefficients)
    set_p(
        image, bandName, "{bandName} * i.cFactor", {"bandName": "i." + bandName}
    )


def brdf(image, bandName, kvolBand, kgeoBand, coefficients):
    """
    Computes the Bidirectional Reflectance Distribution Function (BRDF) for an Earth Engine
    image and adds it as a new band.

    Args:
        image: The Earth Engine image to be processed.
        bandName: The name of the new band that will store the BRDF values.
        kvolBand: The name of the band representing volumetric kernel.
        kgeoBand: The name of the band representing geometric kernel.
        coefficients: A dictionary containing BRDF coefficients (fiso, fvol, fgeo).

    Examples:
        # Compute and add BRDF for a band
        coefficients = {'fiso': 0.0774, 'fgeo': 0.0079, 'fvol': 0.0372}
        brdf(eeImage, 'brdfBand', 'kvolBand', 'kgeoBand', coefficients)
    """
    # Merge coefficients with kernel band names
    args = merge(
        coefficients,
        {
            "kvol": "3 * i." + kvolBand,
            "kgeo": "i." + kgeoBand,
        },
    )
    # Compute and add the BRDF band
    return set_p(
        image, bandName, "{fiso} + {fvol} * {kvol} + {fgeo} * {kvol}", args
    )


def sunZenOut(image):
    """
    Calculates the sun zenith angle for the center of an Earth Engine image footprint.

    Args:
        image: The Earth Engine image to be processed.

    Examples:
        # Calculate sun zenith angle
        sunZenOut(eeImage)
    """
    # Calculate center latitude in radians and compute sun zenith angle
    set_p(
        image,
        "centerLat",
        ee.Number(
            ee.Geometry(image.get("system:footprint"))
            .bounds()
            .centroid(30)
            .coordinates()
            .get(0)
        )
        .multiply(math.pi)
        .divide(180),
    )
    set_p(
        image,
        "sunZenOut",
        "(31.0076"
        + "- 0.1272 * i.centerLat"
        + "+ 0.01187 * pow(i.centerLat, 2)"
        + "+ 2.40E-05 * pow(i.centerLat, 3)"
        + "- 9.48E-07 * pow(i.centerLat, 4)"
        + "- 1.95E-09 * pow(i.centerLat, 5)"
        + "+ 6.15E-11 * pow(i.centerLat, 6)) * {pi}/180",
    )


def rossThick(image, bandName, sunZen, viewZen, relativeSunViewAz):
    """
    Applies the Ross Thick kernel to a specified band in an Earth Engine image.

    Args:
        image: The Earth Engine image to be processed.
        bandName: The name of the band to which the kernel is applied.
        sunZen: The name of the band representing the solar zenith angle.
        viewZen: The name of the band representing the view zenith angle.
        relativeSunViewAz: The name of the band representing the relative azimuth angle.

    Examples:
        # Apply Ross Thick kernel to a band
        rossThick(eeImage, 'rossThickBand', 'sunZen', 'viewZen', 'relAzimuth')
    """
    # Define arguments for the kernel computation
    args = {
        "sunZen": sunZen,
        "viewZen": viewZen,
        "relativeSunViewAz": relativeSunViewAz,
    }
    # Compute phase angle and Ross Thick kernel
    cosPhaseAngle(image, "cosPhaseAngle", sunZen, viewZen, relativeSunViewAz)
    set_p(image, "phaseAngle", "acos(i.cosPhaseAngle)")
    set_p(
        image,
        bandName,
        "(({pi}/2 - i.phaseAngle) * i.cosPhaseAngle + sin(i.phaseAngle)) "
        + "/ (cos({sunZen}) + cos({viewZen})) - {pi}/4",
        args,
    )


def solarPosition(image):
    """
    Calculates the solar position parameters for an Earth Engine image.

    Args:
        image: The Earth Engine image to be processed.

    Examples:
        # Calculate solar position parameters
        solarPosition(eeImage)
    """
    # Extract and process date and time
    date = ee.Date(ee.Number(image.get("system:time_start")))
    secondsInHour = 3600
    set_p(image, "longDeg", ee.Image.pixelLonLat().select("longitude"))
    set_p(
        image,
        "latRad",
        ee.Image.pixelLonLat().select("latitude").multiply(math.pi).divide(180),
    )
    set_p(
        image,
        "hourGMT",
        ee.Number(date.getRelative("second", "day")).divide(secondsInHour),
    )
    # Calculate Julian Date Proportion and related parameters
    set_p(image, "jdp", date.getFraction("year"))  # Julian Date Proportion
    set_p(
        image, "jdpr", "i.jdp * 2 * {pi}"  # Julian Date Proportion in Radians
    )
    set_p(image, "meanSolarTime", "i.hourGMT + i.longDeg / 15")
    set_p(
        image,
        "localSolarDiff",
        "(0.000075 + 0.001868 * cos(i.jdpr) - 0.032077 * sin(i.jdpr)"
        + "- 0.014615 * cos(2 * i.jdpr) - 0.040849 * sin(2 * i.jdpr))"
        + "* 12 * 60 / {pi}",
    )
    set_p(
        image, "trueSolarTime", "i.meanSolarTime + i.localSolarDiff / 60 - 12"
    )
    set_p(image, "angleHour", "i.trueSolarTime * 15 * {pi} / 180")
    set_p(
        image,
        "delta",
        "0.006918 - 0.399912 * cos(i.jdpr) + 0.070257 * sin(i.jdpr) - 0.006758 * cos(2 * i.jdpr)"
        + "+ 0.000907 * sin(2 * i.jdpr) - 0.002697 * cos(3 * i.jdpr) + 0.001480 * sin(3 * i.jdpr)",
    )
    set_p(
        image,
        "cosSunZen",
        "sin(i.latRad) * sin(i.delta) "
        + "+ cos(i.latRad) * cos(i.delta) * cos(i.angleHour)",
    )
    set_p(image, "sunZen", "acos(i.cosSunZen)")
    set_p(
        image,
        "sinSunAzSW",
        to_image("cos(i.delta) * sin(i.angleHour) / sin(i.sunZen)").clamp(
            -1, 1
        ),
    )
    set_p(
        image,
        "cosSunAzSW",
        "(-cos(i.latRad) * sin(i.delta)"
        + "+ sin(i.latRad) * cos(i.delta) * cos(i.angleHour)) / sin(i.sunZen)",
    )
    set_p(image, "sunAzSW", "asin(i.sinSunAzSW)")
    set_pIf(
        image, "sunAzSW", "i.cosSunAzSW <= 0", "{pi} - i.sunAzSW", "sunAzSW"
    )
    set_pIf(
        image,
        "sunAzSW",
        "i.cosSunAzSW > 0 and i.sinSunAzSW <= 0",
        "2 * {pi} + i.sunAzSW",
        "sunAzSW",
    )
    set_p(image, "sunAz", "i.sunAzSW + {pi}")
    set_pIf(image, "sunAz", "i.sunAz > 2 * {pi}", "i.sunAz - 2 * {pi}", "sunAz")


def viewAngles(image, corners):
    """
    Calculates the view angles for an Earth Engine image based on its corner coordinates.

    Args:
        image: The Earth Engine image to be processed.
        corners: A dictionary containing the corner coordinates of the image.

    Examples:
        # Calculate view angles
        corners = findCorners(eeImage)
        viewAngles(eeImage, corners)
    """
    maxDistanceToSceneEdge = 1000000
    maxSatelliteZenith = 7.5
    # Calculate upper and lower center points
    upperCenter = pointBetween(
        corners.get("upperLeft"), corners.get("upperRight")
    )
    lowerCenter = pointBetween(
        corners.get("lowerLeft"), corners.get("lowerRight")
    )
    # Calculate slope and view azimuth
    slope = slopeBetween(lowerCenter, upperCenter)
    slopePerp = ee.Number(-1).divide(slope)
    set_p(
        image,
        "viewAz",
        ee.Image(ee.Number(math.pi / 2).subtract((slopePerp).atan())),
    )

    # Calculate view zenith
    leftLine = toLine(corners.get("upperLeft"), corners.get("lowerLeft"))
    rightLine = toLine(corners.get("upperRight"), corners.get("lowerRight"))
    leftDistance = ee.FeatureCollection(leftLine).distance(
        maxDistanceToSceneEdge
    )
    rightDistance = ee.FeatureCollection(rightLine).distance(
        maxDistanceToSceneEdge
    )
    viewZenith = (
        rightDistance.multiply(maxSatelliteZenith * 2)
        .divide(rightDistance.add(leftDistance))
        .subtract(maxSatelliteZenith)
    )
    set_p(image, "viewZen", viewZenith.multiply(math.pi).divide(180))
    
    
def apply_scale_factors(image):
    optical_bands = image.select('SR_B.*') \
        .multiply(0.0000275) \
        .add(-0.2) \
        .multiply(10000) \
        .toInt16()
    return image.addBands(optical_bands, None, True)

    
if __name__ == "__main__":
    
    # Load in a given image
    image = ee.Image('LANDSAT/LC08/C02/T1_L2/LC08_044034_20140318')

    # Apply the scaling factors
    image = apply_scale_factors(image)


