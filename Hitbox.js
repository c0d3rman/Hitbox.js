EPSILON = 0.00001;

class Hitbox {
  constructor(draw = true) {
    this.draw = draw;
    this.components = [];
    this.collisionFnMap = {};

    this.registerCollisionFn("coord", "coord", Hitbox.collideCoordCoord);
    this.registerCollisionFn("coord", "line", Hitbox.collideCoordLine);
    this.registerCollisionFn("coord", "poly", Hitbox.collideCoordPoly);
    this.registerCollisionFn("coord", "ellipse", Hitbox.collideCoordEllipse);
    this.registerCollisionFn("line", "line", Hitbox.collideLineLine);
    this.registerCollisionFn("line", "poly", Hitbox.collideLinePoly);
    this.registerCollisionFn("line", "ellipse", Hitbox.collideLineEllipse);
    this.registerCollisionFn("poly", "poly", Hitbox.collidePolyPoly);
    this.registerCollisionFn("poly", "ellipse", Hitbox.collidePolyEllipse);
    this.registerCollisionFn(
      "ellipse",
      "ellipse",
      Hitbox.collideEllipseEllipse
    );
  }

  // Public functions - collisions
  // =============================

  doesCollide(other) {
    for (let i = 0; i < this.components.length; i++) {
      for (let j = 0; j < other.components.length; j++) {
        let key = this.components[i].type + "," + other.components[j].type;
        let collisionFn = this.collisionFnMap[key];
        if (collisionFn(this.components[i], other.components[j])) {
          return true;
        }
      }
    }
    return false;
  }

  doesCollidePoint(x, y) {
    return this.doesCollide({ components: [{ type: "coord", x: x, y: y }] });
  }

  // Public functions - p5 shape duplications
  // ========================================

  // Special - does not draw
  coord(x, y) {
    [x, y] = Hitbox.transformPoint(x, y);
    this.components.push({
      type: "coord",
      x: x,
      y: y,
    });
  }

  line(x1, y1, x2, y2) {
    if (this.draw) {
      line(x1, y1, x2, y2);
    }
    [x1, y1] = Hitbox.transformPoint(x1, y1);
    [x2, y2] = Hitbox.transformPoint(x2, y2);
    this.components.push(Hitbox.makeLine(x1, y1, x2, y2));
  }

  rect(x, y, width, height) {
    if (this.draw) {
      rect(x, y, width, height);
    }
    this.components.push(
      Hitbox.transformPolygon({
        type: "poly",
        x: [x, x + width, x + width, x],
        y: [y, y, y + height, y + height],
      })
    );
  }

  triangle(x1, y1, x2, y2, x3, y3) {
    if (this.draw) {
      triangle(x1, y1, x2, y2, x3, y3);
    }
    this.components.push(
      Hitbox.transformPolygon({
        type: "poly",
        x: [x1, x2, x3],
        y: [y1, y2, y3],
      })
    );
  }

  square(x, y, s) {
    this.rect(x, y, s, s);
  }

  quad(x1, y1, x2, y2, x3, y3, x4, y4) {
    if (this.draw) {
      quad(x1, y1, x2, y2, x3, y3, x4, y4);
    }
    this.components.push(
      Hitbox.transformPolygon({
        type: "poly",
        x: [x1, x2, x3, x4],
        y: [y1, y2, y3, y4],
      })
    );
  }

  ellipse(x, y, width, height) {
    if (this.draw) {
      ellipse(x, y, width, height);
    }

    // Get conjugate diameters by transforming the original vertices of the ellipse
    let [vertex1x, vertex1y] = Hitbox.transformPoint(x + width / 2, y);
    let [vertex2x, vertex2y] = Hitbox.transformPoint(x, y + height / 2);
    let [centerX, centerY] = Hitbox.transformPoint(x, y);

    this.components.push(
      Hitbox.makeEllipse(
        centerX,
        centerY,
        vertex1x,
        vertex1y,
        vertex2x,
        vertex2y
      )
    );
  }

  circle(x, y, d) {
    this.ellipse(x, y, d, d);
  }

  // Private helpers
  // ===============

  registerCollisionFn(type1, type2, fn) {
    this.collisionFnMap[type1 + "," + type2] = fn;
    if (type1 != type2) {
      this.collisionFnMap[type2 + "," + type1] = function (a, b) {
        return fn(b, a);
      };
    }
  }

  static transformPoint(x, y, mat = null) {
    let scaleFactor = 1;
    if (mat == null) {
      mat = drawingContext.getTransform();
      scaleFactor = pixelDensity();
    }

    let point = mat.transformPoint(new DOMPoint(x, y));
    return [point.x / scaleFactor, point.y / scaleFactor];
  }

  static transformLine(line, mat = null) {
    let [x1Rot, y1Rot] = Hitbox.transformPoint(line.x1, line.y1, mat);
    let [x2Rot, y2Rot] = Hitbox.transformPoint(line.x2, line.y2, mat);
    return Hitbox.makeLine(x1Rot, y1Rot, x2Rot, y2Rot);
  }

  static transformPolygon(poly, mat = null) {
    let tpoly = Object.assign({}, poly);
    for (let i = 0; i < poly.x.length; i++) {
      let [xRot, yRot] = Hitbox.transformPoint(poly.x[i], poly.y[i], mat);
      tpoly.x[i] = xRot;
      tpoly.y[i] = yRot;
    }
    return tpoly;
  }

  static transformEllipse(ellipse, mat = null) {
    let [xRot, yRot] = Hitbox.transformPoint(ellipse.x, ellipse.y, mat);
    let [dia1xRot, dia1yRot] = Hitbox.transformPoint(
      ellipse.dia1x,
      ellipse.dia1y,
      mat
    );
    let [dia2xRot, dia2yRot] = Hitbox.transformPoint(
      ellipse.dia2x,
      ellipse.dia2y,
      mat
    );
    return Hitbox.makeEllipse(
      xRot,
      yRot,
      dia1xRot,
      dia1yRot,
      dia2xRot,
      dia2yRot
    );
  }

  static distance(x1, y1, x2, y2) {
    return Math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2);
  }

  static makeLine(x1, y1, x2, y2) {
    let mag = Hitbox.distance(x1, y1, x2, y2);
    let vx = x2 - x1;
    let vy = y2 - y1;
    if (mag != 0) {
      vx /= mag;
      vy /= mag;
    }
    let t;
    if (vx != 0) {
      t = (x2 - x1) / vx;
    } else {
      t = (y2 - y1) / vy;
    }
    return {
      type: "line",
      x1: x1,
      y1: y1,
      x2: x2,
      y2: y2,
      vx: vx,
      vy: vy,
      t: t,
    };
  }

  // Makes an ellipse from two conjugate diameters, defined by a center and two points on the ellipse
  // This math is from wikipedia: https://en.wikipedia.org/wiki/Ellipse#General_ellipse_2
  static makeEllipse(centerX, centerY, dia1x, dia1y, dia2x, dia2y) {
    let vf1x = dia1x - centerX;
    let vf1y = dia1y - centerY;
    let vf2x = dia2x - centerX;
    let vf2y = dia2y - centerY;

    // Get the t of the first vertex (second vertex is just this plus PI/2)
    let denom = 2 * (vf1x * vf2x + vf1y * vf2y);
    let t0;
    if (denom != 0) {
      t0 =
        (Math.PI / 2 -
          Math.atan(
            (vf1x * vf1x + vf1y * vf1y - (vf2x * vf2x + vf2y * vf2y)) / denom
          )) /
        2;
    } else {
      // If the diameters are perpendicular, then t0 = 0
      t0 = 0;
    }

    let v1x = centerX + vf1x * Math.cos(t0) + vf2x * Math.sin(t0);
    let v1y = centerY + vf1y * Math.cos(t0) + vf2y * Math.sin(t0);
    let v2x =
      centerX + vf1x * Math.cos(t0 + PI / 2) + vf2x * Math.sin(t0 + PI / 2);
    let v2y =
      centerY + vf1y * Math.cos(t0 + PI / 2) + vf2y * Math.sin(t0 + PI / 2);

    // Calculate semiaxes using vertices
    let a = Hitbox.distance(centerX, centerY, v1x, v1y);
    let b = Hitbox.distance(centerX, centerY, v2x, v2y);

    // Make sure a is the major and b is the minor
    if (b > a) {
      [a, b] = [b, a];
      [v1x, v1y, v2x, v2y] = [v2x, v2y, v1x, v1y];
    }

    // Calculate angle of rotation
    let angle = Math.atan2(v1y - centerY, v1x - centerX);

    // Calculate foci using semiaxes
    let c = Math.sqrt(a ** 2 - b ** 2);
    let f1x = centerX + c * Math.cos(angle);
    let f1y = centerY + c * Math.sin(angle);
    let f2x = centerX - c * Math.cos(angle);
    let f2y = centerY - c * Math.sin(angle);

    return {
      type: "ellipse",
      x: centerX,
      y: centerY,
      dia1x: dia1x,
      dia1y: dia1y,
      dia2x: dia2x,
      dia2y: dia2y,
      a: a,
      b: b,
      f1x: f1x,
      f1y: f1y,
      f2x: f2x,
      f2y: f2y,
      r: a * 2,
      angle: angle,
    };
  }

  // Check which side of a line a point is on
  static sideOfLine(x, y, line) {
    return Math.sign(
      line.vy * x - line.vx * y - (line.vy * line.x1 - line.vx * line.y1)
    );
  }

  // Private static collision functions
  // ==================================

  static collideCoordCoord(coord1, coord2) {
    return (
      Math.abs(coord1.x - coord2.x) < EPSILON &&
      Math.abs(coord1.y - coord2.y) < EPSILON
    );
  }

  static collideCoordLine(coord, line) {
    return (
      Math.abs(
        Hitbox.distance(coord.x, coord.y, line.x1, line.y1) +
          Hitbox.distance(coord.x, coord.y, line.x2, line.y2) -
          Hitbox.distance(line.x1, line.y1, line.x2, line.y2)
      ) < EPSILON
    );
  }

  static collideCoordPoly(coord, poly) {
    // Ray casting algorithm - see https://stackoverflow.com/a/218081/2674563

    // Pick an arbitrary point outside the polygon
    let x = min(poly.x) - 1000;
    let y = min(poly.y) - 1000;

    // Draw a line between the point and the target coord
    let line = Hitbox.makeLine(x, y, coord.x, coord.y);

    // Count how many times the line intersects the polygon

    let intersections = 0;
    for (let i = 0; i < poly.x.length; i++) {
      if (
        Hitbox.collideLineLine(
          line,
          Hitbox.makeLine(
            poly.x[i],
            poly.y[i],
            poly.x[(i + 1) % poly.x.length],
            poly.y[(i + 1) % poly.x.length]
          )
        )
      ) {
        intersections++;
      }

      // If we intersect with any vertices, then we double-count, so subtract back
      // Also need to check sidedness - see https://stackoverflow.com/q/14130742/2674563
      if (
        Hitbox.sideOfLine(
          poly.x[(i + 1) % poly.x.length],
          poly.y[(i + 1) % poly.x.length],
          line
        ) !=
          Hitbox.sideOfLine(
            poly.x[(i - 1 + poly.x.length) % poly.x.length],
            poly.y[(i - 1 + poly.x.length) % poly.x.length],
            line
          ) &&
        Hitbox.collideCoordLine({ x: poly.x[i], y: poly.y[i] }, line)
      ) {
        intersections--;
      }
    }

    // If the number of intersections is odd, the target coord is inside the polygon
    return intersections % 2 == 1;
  }

  static collideCoordEllipse(coord, ellipse) {
    return (
      Hitbox.distance(coord.x, coord.y, ellipse.f1x, ellipse.f1y) +
        Hitbox.distance(coord.x, coord.y, ellipse.f2x, ellipse.f2y) <=
      ellipse.r
    );
  }

  static collideLineLine(line1, line2) {
    let denom = line1.vx * line2.vy - line2.vx * line1.vy;

    if (denom == 0) {
      // Lines are parallel
      return (
        Hitbox.collideCoordLine({ x: line1.x1, y: line1.y1 }, line2) ||
        Hitbox.collideCoordLine({ x: line1.x2, y: line1.y2 }, line2)
      );
    }

    let t1 =
      (line2.vx * (line1.y1 - line2.y1) + line2.vy * (line2.x1 - line1.x1)) /
      denom;
    let intersectionX = line1.x1 + line1.vx * t1;
    let intersectionY = line1.y1 + line1.vy * t1;

    return (
      Hitbox.collideCoordLine({ x: intersectionX, y: intersectionY }, line1) &&
      Hitbox.collideCoordLine({ x: intersectionX, y: intersectionY }, line2)
    );
  }

  static collideLinePoly(line, poly) {
    return (
      Hitbox.collideCoordPoly({ x: line.x1, y: line.y1 }, poly) ||
      Hitbox.collideCoordPoly({ x: line.x2, y: line.y2 }, poly) ||
      [...Array(poly.x.length).keys()].reduce(
        (l, i) =>
          l ||
          Hitbox.collideLineLine(
            Hitbox.makeLine(
              poly.x[i],
              poly.y[i],
              poly.x[(i + 1) % poly.x.length],
              poly.y[(i + 1) % poly.x.length]
            ),
            line
          )
      )
    );
  }

  static collideLineEllipse(line, ellipse) {
    // Translate and rotate the line so we can put the ellipse flat at the origin, for cleaner math
    let mat = new DOMMatrix()
      .rotate((180 / PI) * -ellipse.angle)
      .translate(-ellipse.x, -ellipse.y);
    let tline = Hitbox.transformLine(line, mat);

    // Calculate discriminant to see if the line intersects the ellipse
    // (Math is from sagemath, gotten by plugging x = x0 + vx * t and y = y0 + vy * t into x^2 / a^2 + y^2 / b^2 = 1 and expanding into a quadratic of t)
    let a = tline.x1 ** 2 / ellipse.a ** 2 + tline.y1 ** 2 / ellipse.b ** 2 - 1;
    let b =
      (2 * tline.vx * tline.x1) / ellipse.a ** 2 +
      (2 * tline.vy * tline.y1) / ellipse.b ** 2;
    let c = tline.vx ** 2 / ellipse.a ** 2 + tline.vy ** 2 / ellipse.b ** 2;
    let disc = b ** 2 - 4 * a * c;

    if (Math.sign(disc) >= 0) {
      // If an intersection with the line occurs, we calculate the t1 and t2 where it happens to see if it lies on our line *segment*
      let t1 =
        -(
          ellipse.b ** 2 * tline.vx * tline.x1 +
          ellipse.a ** 2 * tline.vy * tline.y1 +
          Math.sqrt(
            ellipse.b ** 2 * tline.vx ** 2 +
              ellipse.a ** 2 * tline.vy ** 2 -
              tline.vy ** 2 * tline.x1 ** 2 +
              2 * tline.vx * tline.vy * tline.x1 * tline.y1 -
              tline.vx ** 2 * tline.y1 ** 2
          ) *
            ellipse.a *
            ellipse.b
        ) /
        (ellipse.b ** 2 * tline.vx ** 2 + ellipse.a ** 2 * tline.vy ** 2);
      let t2 =
        -(
          ellipse.b ** 2 * tline.vx * tline.x1 +
          ellipse.a ** 2 * tline.vy * tline.y1 -
          Math.sqrt(
            ellipse.b ** 2 * tline.vx ** 2 +
              ellipse.a ** 2 * tline.vy ** 2 -
              tline.vy ** 2 * tline.x1 ** 2 +
              2 * tline.vx * tline.vy * tline.x1 * tline.y1 -
              tline.vx ** 2 * tline.y1 ** 2
          ) *
            ellipse.a *
            ellipse.b
        ) /
        (ellipse.b ** 2 * tline.vx ** 2 + ellipse.a ** 2 * tline.vy ** 2);
      if ((t1 >= 0 && t1 <= tline.t) || (t2 >= 0 && t2 <= tline.t)) {
        return true;
      }
    }

    // If the line doesn't intersect the ellipse, then we just check its endpoints - it could still be entirely within the ellipse
    // Technically this could equivalently be && instead of || but no reason not to play it safe
    return (
      Hitbox.collideCoordEllipse({ x: line.x1, y: line.y1 }, ellipse) ||
      Hitbox.collideCoordEllipse({ x: line.x2, y: line.y2 }, ellipse)
    );
  }

  static collidePolyPoly(poly1, poly2) {
    for (let i = 0; i < poly1.x.length; i++) {
      let line1 = Hitbox.makeLine(
        poly1.x[i],
        poly1.y[i],
        poly1.x[(i + 1) % poly1.x.length],
        poly1.y[(i + 1) % poly1.x.length]
      );
      for (let j = 0; j < poly2.x.length; j++) {
        let line2 = Hitbox.makeLine(
          poly2.x[j],
          poly2.y[j],
          poly2.x[(j + 1) % poly2.x.length],
          poly2.y[(j + 1) % poly2.x.length]
        );
        if (Hitbox.collideLineLine(line1, line2)) {
          return true;
        }
      }
    }

    return (
      [...Array(poly1.x.length).keys()].reduce(
        (l, i) =>
          l || Hitbox.collideCoordPoly({ x: poly1.x[i], y: poly1.y[i] }, poly2)
      ) ||
      [...Array(poly2.x.length).keys()].reduce(
        (l, i) =>
          l || Hitbox.collideCoordPoly({ x: poly2.x[i], y: poly2.y[i] }, poly1)
      )
    );
  }

  static collidePolyEllipse(poly, ellipse) {
    for (let i = 0; i < poly.x.length; i++) {
      let line = Hitbox.makeLine(
        poly.x[i],
        poly.y[i],
        poly.x[(i + 1) % poly.x.length],
        poly.y[(i + 1) % poly.x.length]
      );

      if (Hitbox.collideLineEllipse(line, ellipse)) {
        return true;
      }
    }

    return (
      [...Array(poly.x.length).keys()].reduce(
        (l, i) =>
          l ||
          Hitbox.collideCoordEllipse({ x: poly.x[i], y: poly.y[i] }, ellipse)
      ) || Hitbox.collideCoordPoly({ x: ellipse.x, y: ellipse.y }, poly)
    );
  }

  // This code is adapted from https://www.khanacademy.org/computer-programming/c/5567955982876672
  static collideEllipseEllipse(ellipse1, ellipse2) {
    /*
     * Does the quartic function described by
     * y = z4*x⁴ + z3*x³ + z2*x² + z1*x + z0 have *any*
     * real solutions? See
     * http://en.wikipedia.org/wiki/Quartic_function
     * Thanks to Dr. David Goldberg for the conversion to
     * a depressed quartic!
     */
    let realRoot = function (z4, z3, z2, z1, z0) {
      /* First trivial checks for z0 or z4 being zero */
      if (Math.abs(z0) < EPSILON) {
        return true; /* zero is a root! */
      }
      if (Math.abs(z4) < EPSILON) {
        if (Math.abs(z3) >= EPSILON) {
          return true; /* cubics always have roots */
        }
        if (Math.abs(z2) >= EPSILON) {
          return z1 ** 2 - 4 * z2 * z0 >= 0; /* quadratic, use determinant */
        }
        return Math.abs(z1) >= EPSILON; /* sloped lines have one root */
      }
      let a = z3 / z4,
        b = z2 / z4,
        c = z1 / z4,
        d = z0 / z4;
      let p = (8 * b - 3 * a ** 2) / 8;
      let q = (a ** 3 - 4 * a * b + 8 * c) / 8;
      let r = (-3 * a ** 4 + 256 * d - 64 * c * a + 16 * a ** 2 * b) / 256;
      /*
       *   x⁴ +        p*x² + q*x + r
       * a*x⁴ + b*x³ + c*x² + d*x + e
       * so a=1  b=0  c=p  d=q  e=r
       * That is, we have a depessed quartic.
       */
      let discrim =
        256 * r ** 3 -
        128 * p ** 2 * r ** 2 +
        144 * p * q ** 2 * r -
        27 * q ** 4 +
        16 * p ** 4 * r -
        4 * p ** 3 * q ** 2;
      let P = 8 * p;
      let D = 64 * r - 16 * p ** 2;

      return (
        discrim < 0 ||
        (discrim > 0 && P < 0 && D < 0) ||
        (discrim === 0 && (D !== 0 || P <= 0))
      );
    };

    /*
     * Is the Y coordinate(s) of the intersection of two conic
     * sections real? They are in their bivariate form,
     * ax²  + bxy  + cx²  + dx  + ey  + f = 0
     * For now, a and a1 cannot be zero.
     * See https://docs.google.com/file/d/0B7wsEy6bpVePSEt2Ql9hY0hFdjA/
     */
    let yIntersect = function (a, b, c, d, e, f, a1, b1, c1, d1, e1, f1) {
      /*
       * Normalize the conics by their first coefficient, a.
       * Then get the differnce of the two equations.
       */
      let deltaB = (b1 /= a1) - (b /= a),
        deltaC = (c1 /= a1) - (c /= a),
        deltaD = (d1 /= a1) - (d /= a),
        deltaE = (e1 /= a1) - (e /= a),
        deltaF = (f1 /= a1) - (f /= a);

      /* Special case for b's and d's being equal */
      if (deltaB === 0 && deltaD === 0) {
        return realRoot(0, 0, deltaC, deltaE, deltaF);
      }

      let a3 = b * c1 - b1 * c;
      let a2 = b * e1 + d * c1 - b1 * e - d1 * c;
      a1 = b * f1 + d * e1 - b1 * f - d1 * e;
      let a0 = d * f1 - d1 * f;

      let A = deltaC * deltaC - a3 * deltaB;
      let B = 2 * deltaC * deltaE - deltaB * a2 - deltaD * a3;
      let C = deltaE * deltaE + 2 * deltaC * deltaF - deltaB * a1 - deltaD * a2;
      let D = 2 * deltaE * deltaF - deltaD * a1 - deltaB * a0;
      let E = deltaF * deltaF - deltaD * a0;
      return realRoot(A, B, C, D, E);
    };

    /*
     * Do two conics sections el and el1 intersect? Each are in
     * bivariate form, ax²  + bxy  + cx²  + dx  + ey  + f = 0
     * Solve by constructing a quartic that must have a real
     * solution if they intersect.  This checks for real Y
     * intersects, then flips the parameters around to check
     * for real X intersects.
     */
    var conicsIntersect = function (el, el1) {
      /* check for real y intersects, then real x intersects */
      return (
        yIntersect(
          el.a,
          el.b,
          el.c,
          el.d,
          el.e,
          el.f,
          el1.a,
          el1.b,
          el1.c,
          el1.d,
          el1.e,
          el1.f
        ) &&
        yIntersect(
          el.c,
          el.b,
          el.a,
          el.e,
          el.d,
          el.f,
          el1.c,
          el1.b,
          el1.a,
          el1.e,
          el1.d,
          el1.f
        )
      );
    };

    /*
     * Express an ellipse in terms of a "bivariate"
     * polynomial that sums to zero. See
     * http://elliotnoma.wordpress.com/2013/04/10/a-closed-form-solution-for-the-intersections-of-two-ellipses
     */
    let bivariateForm = function (ellipse) {
      let A = Math.cos(ellipse.angle);
      let B = Math.sin(ellipse.angle);
      /*
       * Start by rotating the ellipse center by the OPPOSITE
       * of the desired angle.  That way when the bivariate
       * computation transforms it back, it WILL be at the
       * correct (and original) coordinates.
       */
      let r = new DOMMatrix()
        .rotate((180 / PI) * -ellipse.angle)
        .transformPoint(new DOMPoint(ellipse.x, ellipse.y));
      let a = r.x,
        c = r.y;

      /*
       * Now let the bivariate computation rotate in the opposite direction.
       * We only need to flip B to do this, because cos(-theta) = cos(theta)
       */
      B = -B;
      let b = ellipse.a ** 2,
        d = ellipse.b ** 2;
      return {
        a: (A * A) / b + (B * B) / d /* x² coefficient */,
        b: (-2 * A * B) / b + (2 * A * B) / d /* xy coeff */,
        c: (B * B) / b + (A * A) / d /* y² coeff */,
        d: (-2 * a * A) / b - (2 * c * B) / d /* x coeff */,
        e: (2 * a * B) / b - (2 * c * A) / d /* y coeff */,
        f: (a * a) / b + (c * c) / d - 1 /* constant */,
        /* So, ax² + bxy + cy² + dx + ey + f = 0 */
      };
    };

    // First check if one ellipse is within the other
    if (
      Hitbox.collideCoordEllipse({ x: ellipse1.x, y: ellipse1.y }, ellipse2) ||
      Hitbox.collideCoordEllipse({ x: ellipse2.x, y: ellipse2.y }, ellipse1)
    ) {
      return true;
    }

    // Convert both ellipses to bivariate form, then check if they intersect
    let elps1 = bivariateForm(ellipse1);
    let elps2 = bivariateForm(ellipse2);
    return conicsIntersect(elps1, elps2);
  }
}
