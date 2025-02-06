exports.apply_brdf = function (image) {
  
    var inputBandNames = image.bandNames();
    var constants = {
      pi: Math.PI
    };
  
    var coefficientsByBand = {
      'B1': {fiso: 0.0774, fgeo: 0.0079, fvol: 0.0372},
      'B2': {fiso: 0.1306, fgeo: 0.0178, fvol: 0.0580},
      'B3': {fiso: 0.1690, fgeo: 0.0227, fvol: 0.0574},
      'B4': {fiso: 0.3093, fgeo: 0.0330, fvol: 0.1535},
      'B5': {fiso: 0.3430, fgeo: 0.0453, fvol: 0.1154},
      'B7': {fiso: 0.2658, fgeo: 0.0387, fvol: 0.0639}
    };
    
    var corners = findCorners();
    
    viewAngles();
    solarPosition();
    sunZenOut();
    set('relativeSunViewAz', 'i.sunAz - i.viewAz');
    rossThick('kvol', 'i.sunZen', 'i.viewZen', 'i.relativeSunViewAz');
    rossThick('kvol0', 'i.sunZenOut', 0, 0);
    liThin('kgeo', 'i.sunZen', 'i.viewZen', 'i.relativeSunViewAz');
    liThin('kgeo0', 'i.sunZenOut', 0, 0);
    adjustBands();
    return image.select(inputBandNames);
    
  
    function viewAngles() {
      var maxDistanceToSceneEdge = 1000000;
      var maxSatelliteZenith = 7.5;
      var upperCenter = pointBetween(corners.upperLeft, corners.upperRight);
      var lowerCenter = pointBetween(corners.lowerLeft, corners.lowerRight);
      var slope = slopeBetween(lowerCenter, upperCenter);
      var slopePerp = ee.Number(-1).divide(slope);
      set('viewAz', ee.Image(ee.Number(Math.PI / 2).subtract((slopePerp).atan())));
  
      var leftLine = toLine(corners.upperLeft, corners.lowerLeft)
      var rightLine = toLine(corners.upperRight, corners.lowerRight)
      var leftDistance = ee.FeatureCollection(leftLine).distance(maxDistanceToSceneEdge)
      var rightDistance = ee.FeatureCollection(rightLine).distance(maxDistanceToSceneEdge)
      var viewZenith = rightDistance.multiply(maxSatelliteZenith * 2)
        .divide(rightDistance.add(leftDistance))
        .subtract(maxSatelliteZenith)
      set('viewZen',
        viewZenith.multiply(Math.PI).divide(180))
    }
  
    function solarPosition() {
      // Ported from http://pythonfmask.org/en/latest/_modules/fmask/landsatangles.html
      var date = ee.Date(ee.Number(image.get('system:time_start')))
      var secondsInHour = 3600
      set('longDeg',
        ee.Image.pixelLonLat().select('longitude'))
      set('latRad',
        ee.Image.pixelLonLat().select('latitude')
          .multiply(Math.PI).divide(180))
      set('hourGMT',
        ee.Number(date.getRelative('second', 'day')).divide(secondsInHour))
      set('jdp', // Julian Date Proportion
        date.getFraction('year'))
      set('jdpr', // Julian Date Proportion in Radians
        'i.jdp * 2 * {pi}')
      set('meanSolarTime',
        'i.hourGMT + i.longDeg / 15')
      set('localSolarDiff',
        '(0.000075 + 0.001868 * cos(i.jdpr) - 0.032077 * sin(i.jdpr)' +
        '- 0.014615 * cos(2 * i.jdpr) - 0.040849 * sin(2 * i.jdpr))' +
        '* 12 * 60 / {pi}')
      set('trueSolarTime',
        'i.meanSolarTime + i.localSolarDiff / 60 - 12')
      set('angleHour',
        'i.trueSolarTime * 15 * {pi} / 180')
      set('delta',
        '0.006918 - 0.399912 * cos(i.jdpr) + 0.070257 * sin(i.jdpr) - 0.006758 * cos(2 * i.jdpr)' +
        '+ 0.000907 * sin(2 * i.jdpr) - 0.002697 * cos(3 * i.jdpr) + 0.001480 * sin(3 * i.jdpr)')
      set('cosSunZen',
        'sin(i.latRad) * sin(i.delta) ' +
        '+ cos(i.latRad) * cos(i.delta) * cos(i.angleHour)')
      set('sunZen',
        'acos(i.cosSunZen)')
      set('sinSunAzSW',
        toImage('cos(i.delta) * sin(i.angleHour) / sin(i.sunZen)')
          .clamp(-1, 1))
      set('cosSunAzSW',
        '(-cos(i.latRad) * sin(i.delta)' +
        '+ sin(i.latRad) * cos(i.delta) * cos(i.angleHour)) / sin(i.sunZen)')
      set('sunAzSW',
        'asin(i.sinSunAzSW)')
      setIf('sunAzSW',
        'i.cosSunAzSW <= 0',
        '{pi} - i.sunAzSW',
        'sunAzSW')
      setIf('sunAzSW',
        'i.cosSunAzSW > 0 and i.sinSunAzSW <= 0',
        '2 * {pi} + i.sunAzSW',
        'sunAzSW')
      set('sunAz',
        'i.sunAzSW + {pi}')
      setIf('sunAz',
        'i.sunAz > 2 * {pi}',
        'i.sunAz - 2 * {pi}',
        'sunAz')
    }
  
    function sunZenOut() {
      // https://nex.nasa.gov/nex/static/media/publication/HLS.v1.0.UserGuide.pdf
      set('centerLat',
        ee.Number(
          ee.Geometry(image.get('system:footprint'))
            .bounds().centroid(30).coordinates().get(0))
          .multiply(Math.PI).divide(180))
      set('sunZenOut',
        '(31.0076' +
        '- 0.1272 * i.centerLat' +
        '+ 0.01187 * pow(i.centerLat, 2)' +
        '+ 2.40E-05 * pow(i.centerLat, 3)' +
        '- 9.48E-07 * pow(i.centerLat, 4)' +
        '- 1.95E-09 * pow(i.centerLat, 5)' +
        '+ 6.15E-11 * pow(i.centerLat, 6)) * {pi}/180')
    }
  
    function rossThick(bandName, sunZen, viewZen, relativeSunViewAz) {
      var args = {sunZen: sunZen, viewZen: viewZen, relativeSunViewAz: relativeSunViewAz}
      cosPhaseAngle('cosPhaseAngle', sunZen, viewZen, relativeSunViewAz)
      set('phaseAngle',
        'acos(i.cosPhaseAngle)')
      set(bandName,
        '(({pi}/2 - i.phaseAngle) * i.cosPhaseAngle + sin(i.phaseAngle)) ' +
        '/ (cos({sunZen}) + cos({viewZen})) - {pi}/4', args)
    }
  
    function liThin(bandName, sunZen, viewZen, relativeSunViewAz) {
      // From https://modis.gsfc.nasa.gov/data/atbd/atbd_mod09.pdf
      var args = {
        sunZen: sunZen,
        viewZen: viewZen,
        relativeSunViewAz: relativeSunViewAz,
        'h/b': 2,
      }
  
      anglePrime('sunZenPrime', sunZen)
      anglePrime('viewZenPrime', viewZen)
      cosPhaseAngle('cosPhaseAnglePrime', 'i.sunZenPrime', 'i.viewZenPrime', relativeSunViewAz)
      set('distance',
        'sqrt(pow(tan(i.sunZenPrime), 2) + pow(tan(i.viewZenPrime), 2)' +
        '- 2 * tan(i.sunZenPrime) * tan(i.viewZenPrime) * cos({relativeSunViewAz}))', args)
      set('temp',
        '1/cos(i.sunZenPrime) + 1/cos(i.viewZenPrime)')
      set('cosT',
        toImage('{h/b} * sqrt(pow(i.distance, 2) + pow(tan(i.sunZenPrime) * tan(i.viewZenPrime) * sin({relativeSunViewAz}), 2))' +
          '/ i.temp', args)
          .clamp(-1, 1))
      set('t', 'acos(i.cosT)')
      set('overlap',
        '(1/{pi}) * (i.t - sin(i.t) * i.cosT) * (i.temp)')
      setIf('overlap', 'i.overlap > 0', 0)
      set(bandName,
        'i.overlap - i.temp' +
        '+ (1/2) * (1 + i.cosPhaseAnglePrime) * (1/cos(i.sunZenPrime)) * (1/cos(i.viewZenPrime))')
    }
  
    function anglePrime(name, angle) {
      var args = {'b/r': 1, angle: angle}
      set('tanAnglePrime',
        '{b/r} * tan({angle})', args)
      setIf('tanAnglePrime', 'i.tanAnglePrime < 0', 0)
      set(name,
        'atan(i.tanAnglePrime)')
    }
  
    function cosPhaseAngle(name, sunZen, viewZen, relativeSunViewAz) {
      var args = {
        sunZen: sunZen,
        viewZen: viewZen,
        relativeSunViewAz: relativeSunViewAz
      }
      set(name,
        toImage('cos({sunZen}) * cos({viewZen})' +
          '+ sin({sunZen}) * sin({viewZen}) * cos({relativeSunViewAz})', args)
          .clamp(-1, 1))
    }
  
    function adjustBands() {
      for (var bandName in coefficientsByBand)
        applyCFactor(bandName, coefficientsByBand[bandName])
    }
  
    function applyCFactor(bandName, coefficients) {
      brdf('brdf', 'kvol', 'kgeo', coefficients)
      brdf('brdf0', 'kvol0', 'kgeo0', coefficients)
      set('cFactor',
        'i.brdf0 / i.brdf', coefficients)
      set(bandName,
        '{bandName} * i.cFactor', {bandName: 'i.' + bandName})
    }
  
    function brdf(bandName, kvolBand, kgeoBand, coefficients) {
      var args = merge(coefficients, {
        // kvol: 'i.' + kvolBand,
        kvol: '3 * i.' + kvolBand,     // check this multiplication factor.  Is there an 'optimal' value?  Without a factor here, there is not enough correction.
        kgeo: 'i.' + kgeoBand
      })
      return set(bandName,
        '{fiso} + {fvol} * {kvol} + {fgeo} * {kvol}', args)
    }
  
    function findCorners() {
      var footprint = ee.Geometry(image.get('system:footprint'))
      var bounds = ee.List(footprint.bounds().coordinates().get(0))
      var coords = footprint.coordinates()
  
      var xs = coords.map(function (item) {
        return x(item)
      })
      var ys = coords.map(function (item) {
        return y(item)
      })
  
      function findCorner(targetValue, values) {
        var diff = values.map(function (value) {
          return ee.Number(value).subtract(targetValue).abs()
        })
        var minValue = diff.reduce(ee.Reducer.min())
        var idx = diff.indexOf(minValue)
        return coords.get(idx)
      }
  
      var lowerLeft = findCorner(x(bounds.get(0)), xs)
      var lowerRight = findCorner(y(bounds.get(1)), ys)
      var upperRight = findCorner(x(bounds.get(2)), xs)
      var upperLeft = findCorner(y(bounds.get(3)), ys)
      return {
        upperLeft: upperLeft,
        upperRight: upperRight,
        lowerRight: lowerRight,
        lowerLeft: lowerLeft
      }
    }
  
    function x(point) {
      return ee.Number(ee.List(point).get(0))
    }
  
    function y(point) {
      return ee.Number(ee.List(point).get(1))
    }
  
    function pointBetween(pointA, pointB) {
      return ee.Geometry.LineString([pointA, pointB]).centroid().coordinates()
    }
  
    function slopeBetween(pointA, pointB) {
      return ((y(pointA)).subtract(y(pointB))).divide((x(pointA)).subtract(x(pointB)))
    }
  
    function toLine(pointA, pointB) {
      return ee.Geometry.LineString([pointA, pointB])
    }
  
  // ************** COMMON HELPERS **************
  
    function set(name, toAdd, args) {
      toAdd = toImage(toAdd, args)
      image = image.addBands(toAdd.rename(name), null, true)
    }
  
    function setIf(name, condition, trueValue, falseValue) {
      condition = toImage(condition)
      var trueMasked = toImage(trueValue).mask(toImage(condition))
      var falseMasked = toImage(falseValue).mask(invertMask(condition))
      var value = trueMasked.unmask(falseMasked)
      set(name, value)
      
  
      function invertMask(mask) {
        return mask.multiply(-1).add(1)
      }
    }
  
    function toImage(band, args) {
      if ((typeof band) === 'string') {
        if (band.indexOf('.') > -1 || band.indexOf(' ') > -1 || band.indexOf('{') > -1) {
          band = image.expression(format(band, args), {i: image})
        } else
          band = image.select(band)
      }
      return ee.Image(band)
    }
  
    function format(s, args) {
      if (!args) args = {}
      var allArgs = merge(constants, args)
      var result = s.replace(/{([^{}]*)}/g,
        function (a, b) {
          var replacement = allArgs[b]
          if (replacement == null) {
            print('Undeclared argument: ' + b, 's: ' + s, args)
            return null
          }
          return allArgs[b]
        }
      )
      if (result.indexOf('{') > -1)
        return format(result, args)
      return result
    }
    
    function merge(o1, o2) {
      function addAll(target, toAdd) {
        for (var key in toAdd) target[key] = toAdd[key]
      }
  
      var result = {}
      addAll(result, o1)
      addAll(result, o2)
      return result
    }
  
  }
    