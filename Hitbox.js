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

  // Public functions - collisions & other
  // =====================================

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
    return this.doesCollide({
      components: [{ type: "coord", p: createVector(x, y) }],
    });
  }

  // Turn drawing on and off. Valid values:
  // true - draw as normal (default)
  // false - don't draw
  // "debug" - draw hitboxes as transparent red shapes
  setDraw(draw) {
    this.draw = draw;
  }

  // Public functions - p5 shape duplications
  // ========================================

  point(...args) {
    this.doDraw(() => point(...args));

    let p;
    if (args.length === 1 && args[0] instanceof p5.Vector) {
      p = args[0];
    } else {
      p = createVector(args[0], args[1]);
    }

    this.components.push({
      type: "coord",
      p: Hitbox.transformPoint(p),
    });
  }

  line(x1, y1, x2, y2) {
    if (this.draw) {
      line(x1, y1, x2, y2);
    }
    this.components.push(
      Hitbox.makeLine(
        Hitbox.transformPoint(createVector(x1, y1)),
        Hitbox.transformPoint(createVector(x2, y2))
      )
    );
  }

  rect(...args) {
    this.doDraw(() => rect(...args));

    let [x, y, width, height] = this.handleRectAndEllipseModes(
      args,
      window._renderer._rectMode
    );

    this.components.push(
      Hitbox.transformPolygon({
        type: "poly",
        ps: [
          createVector(x, y),
          createVector(x + width, y),
          createVector(x + width, y + height),
          createVector(x, y + height),
        ],
      })
    );
  }

  triangle(x1, y1, x2, y2, x3, y3) {
    this.doDraw(() => triangle(x1, y1, x2, y2, x3, y3));

    this.components.push(
      Hitbox.transformPolygon({
        type: "poly",
        ps: [createVector(x1, y1), createVector(x2, y2), createVector(x3, y3)],
      })
    );
  }

  square(x, y, s) {
    this.rect(x, y, s, s);
  }

  quad(x1, y1, x2, y2, x3, y3, x4, y4) {
    this.doDraw(() => quad(x1, y1, x2, y2, x3, y3, x4, y4));

    this.components.push(
      Hitbox.transformPolygon({
        type: "poly",
        ps: [
          createVector(x1, y1),
          createVector(x2, y2),
          createVector(x3, y3),
          createVector(x4, y4),
        ],
      })
    );
  }

  ellipse(...args) {
    this.doDraw(() => ellipse(...args));

    let [x, y, width, height] = this.handleRectAndEllipseModes(
      args,
      window._renderer._ellipseMode
    );

    // Get conjugate diameters by transforming the original vertices of the ellipse
    this.components.push(
      Hitbox.makeEllipse(
        Hitbox.transformPoint(createVector(x, y)),
        Hitbox.transformPoint(createVector(x + width / 2, y)),
        Hitbox.transformPoint(createVector(x, y + height / 2))
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

  handleRectAndEllipseModes(args, mode) {
    // Support for no-height syntax
    if (args.length == 3) {
      args.push(args[2]);
    }

    switch (mode) {
      case p5.prototype.CORNER:
        return args.slice(0, 4);
      case p5.prototype.CORNERS:
        return [
          Math.min(args[0], args[2]),
          Math.min(args[1], args[3]),
          Math.abs(args[2] - args[0]),
          Math.abs(args[3] - args[1]),
        ];
      case p5.prototype.RADIUS:
        return [args[0] - args[2], args[1] - args[3], args[2] * 2, args[3] * 2];
      case p5.prototype.CENTER:
        return [args[0] - args[2] / 2, args[1] - args[3] / 2, args[2], args[3]];
    }
  }

  // A helper function for handling drawing
  doDraw(f) {
    if (this.draw) {
      push();
      if (this.draw == "debug") {
        noStroke();
        fill(color(255, 0, 0, 100));
      }
      f();
      pop();
    }
  }

  static transformPoint(p, mat = null) {
    let scaleFactor = 1;
    if (mat == null) {
      mat = drawingContext.getTransform();
      scaleFactor = pixelDensity();
    }

    let pt = mat.transformPoint(new DOMPoint(p.x, p.y));
    return createVector(pt.x / scaleFactor, pt.y / scaleFactor);
  }

  static transformLine(line, mat = null) {
    return Hitbox.makeLine(
      Hitbox.transformPoint(line.p1, mat),
      Hitbox.transformPoint(line.p2, mat)
    );
  }

  static transformPolygon(poly, mat = null) {
    return {
      type: "poly",
      ps: poly.ps.map((p) => Hitbox.transformPoint(p, mat)),
    };
  }

  static transformEllipse(ellipse, mat = null) {
    return Hitbox.makeEllipse(
      Hitbox.transformPoint(ellipse.center, mat),
      Hitbox.transformPoint(ellipse.dia1, mat),
      Hitbox.transformPoint(ellipse.dia2, mat)
    );
  }

  static makeLine(p1, p2) {
    return {
      type: "line",
      p1: p1,
      p2: p2,
      v: p5.Vector.sub(p2, p1).normalize(),
      t: p2.dist(p1),
    };
  }

  // Makes an ellipse from two conjugate diameters, defined by a center and two points on the ellipse
  // This math is from wikipedia: https://en.wikipedia.org/wiki/Ellipse#General_ellipse_2
  static makeEllipse(center, dia1, dia2) {
    let vf1 = p5.Vector.sub(dia1, center);
    let vf2 = p5.Vector.sub(dia2, center);

    // Get the t of the first vertex
    let denom = 2 * vf1.dot(vf2);
    let t0;
    if (denom != 0) {
      t0 = (Math.PI / 2 - Math.atan((vf1.magSq() - vf2.magSq()) / denom)) / 2;
    } else {
      // If the diameters are perpendicular, then t0 = 0
      t0 = 0;
    }
    let t1 = t0 + PI / 2; // second vertex is just this plus PI/2

    // Get the vertices using t0
    let v1 = p5.Vector.add(
      center,
      p5.Vector.add(
        p5.Vector.mult(vf1, Math.cos(t0)),
        p5.Vector.mult(vf2, Math.sin(t0))
      )
    );
    let v2 = p5.Vector.add(
      center,
      p5.Vector.add(
        p5.Vector.mult(vf1, Math.cos(t1)),
        p5.Vector.mult(vf2, Math.sin(t1))
      )
    );

    // Calculate semiaxes using vertices
    let a = center.dist(v1);
    let b = center.dist(v2);

    // Make sure a is the major and b is the minor
    if (b > a) {
      [a, b] = [b, a];
      [v1, v2] = [v2, v1];
    }

    // Calculate angle of rotation
    let angle = Math.atan2(v1.y - center.y, v1.x - center.x);

    // Calculate foci using semiaxes
    let c = Math.sqrt(a ** 2 - b ** 2);
    // Thankfully fromAngle is unaffected by angleMode
    let f1 = p5.Vector.add(center, p5.Vector.fromAngle(angle, c));
    let f2 = p5.Vector.sub(center, p5.Vector.fromAngle(angle, c));

    return {
      type: "ellipse",
      center: center,
      dia1: dia1,
      dia2: dia2,
      a: a,
      b: b,
      f1: f1,
      f2: f2,
      r: a * 2,
      angle: angle,
    };
  }

  // Check which side of a line a point is on
  static sideOfLine(p, line) {
    return Math.sign(
      line.v.y * p.x -
        line.v.x * p.y -
        (line.v.y * line.p1.x - line.v.x * line.p1.y)
    );
  }

  // Private static collision functions
  // ==================================

  static collideCoordCoord(coord1, coord2) {
    return (
      Math.abs(coord1.p.x - coord2.p.x) < EPSILON &&
      Math.abs(coord1.p.y - coord2.p.y) < EPSILON
    );
  }

  static collideCoordLine(coord, line) {
    return (
      Math.abs(
        coord.p.dist(line.p1) + coord.p.dist(line.p2) - line.p1.dist(line.p2)
      ) < EPSILON
    );
  }

  static collideCoordPoly(coord, poly) {
    // Ray casting algorithm - see https://stackoverflow.com/a/218081/2674563

    // Pick an arbitrary point outside the polygon
    let p = createVector(
      min(poly.ps.map((p) => p.x)) - 10000,
      min(poly.ps.map((p) => p.y)) - 10000
    );

    // Draw a line between the point and the target coord
    let line = Hitbox.makeLine(p, coord.p);

    // Count how many times the line intersects the polygon
    let intersections = 0;
    for (let i = 0; i < poly.ps.length; i++) {
      if (
        Hitbox.collideLineLine(
          line,
          Hitbox.makeLine(poly.ps[i], poly.ps[(i + 1) % poly.ps.length])
        )
      ) {
        intersections++;
      }

      // If we intersect with any vertices, then we double-count, so subtract back
      // Also need to check sidedness - see https://stackoverflow.com/q/14130742/2674563
      if (
        Hitbox.sideOfLine(poly.ps[(i + 1) % poly.ps.length], line) !=
          Hitbox.sideOfLine(
            poly.ps[(i - 1 + poly.ps.length) % poly.ps.length],
            line
          ) &&
        Hitbox.collideCoordLine({ type: "coord", p: poly.ps[i] }, line)
      ) {
        intersections--;
      }
    }

    // If the number of intersections is odd, the target coord is inside the polygon
    return intersections % 2 == 1;
  }

  static collideCoordEllipse(coord, ellipse) {
    // Geometric definition of an ellipse:
    // The set of all points for which the sum of a point's distance from the two foci is equal to some constant
    return ellipse.f1.dist(coord.p) + ellipse.f2.dist(coord.p) <= ellipse.r;
  }

  static collideLineLine(line1, line2) {
    let denom = line1.v.x * line2.v.y - line2.v.x * line1.v.y;

    if (denom == 0) {
      // Lines are parallel
      return (
        Hitbox.collideCoordLine({ type: "coord", p: line1.p1 }, line2) ||
        Hitbox.collideCoordLine({ type: "coord", p: line1.p2 }, line2)
      );
    }

    let t1 =
      (line2.v.x * (line1.p1.y - line2.p1.y) +
        line2.v.y * (line2.p1.x - line1.p1.x)) /
      denom;
    let intersection = p5.Vector.add(line1.p1, p5.Vector.mult(line1.v, t1));

    return (
      Hitbox.collideCoordLine({ type: "coord", p: intersection }, line1) &&
      Hitbox.collideCoordLine({ type: "coord", p: intersection }, line2)
    );
  }

  static collideLinePoly(line, poly) {
    return (
      Hitbox.collideCoordPoly({ type: "coord", p: line.p1 }, poly) ||
      Hitbox.collideCoordPoly({ type: "coord", p: line.p2 }, poly) ||
      [...Array(poly.ps.length).keys()].reduce(
        (l, i) =>
          l ||
          Hitbox.collideLineLine(
            Hitbox.makeLine(poly.ps[i], poly.ps[(i + 1) % poly.ps.length]),
            line
          )
      )
    );
  }

  static collideLineEllipse(line, ellipse) {
    // Translate and rotate the line so we can put the ellipse flat at the origin, for cleaner math
    let mat = new DOMMatrix()
      .rotate((180 / PI) * -ellipse.angle)
      .translate(-ellipse.center.x, -ellipse.center.y);
    let tline = Hitbox.transformLine(line, mat);

    // Calculate discriminant to see if the line intersects the ellipse
    // (Math is from sagemath, gotten by plugging x = x0 + vx * t and y = y0 + vy * t into x^2 / a^2 + y^2 / b^2 = 1 and expanding into a quadratic of t)
    let a =
      tline.p1.x ** 2 / ellipse.a ** 2 + tline.p1.y ** 2 / ellipse.b ** 2 - 1;
    let b =
      (2 * tline.v.x * tline.p1.x) / ellipse.a ** 2 +
      (2 * tline.v.y * tline.p1.y) / ellipse.b ** 2;
    let c = tline.v.x ** 2 / ellipse.a ** 2 + tline.v.y ** 2 / ellipse.b ** 2;
    let disc = b ** 2 - 4 * a * c;

    // If discriminant is positive we have 2 solutions, if it's zero we have 1 solution.
    // Only when it's negative are there no solutions and hence no intersections
    if (Math.sign(disc) >= 0) {
      // If an intersection with the line occurs, we calculate the t1 and t2 where it happens to see if it lies on our line *segment* (as opposed to just our line)
      let t1 =
        -(
          ellipse.b ** 2 * tline.v.x * tline.p1.x +
          ellipse.a ** 2 * tline.v.y * tline.p1.y +
          Math.sqrt(
            ellipse.b ** 2 * tline.v.x ** 2 +
              ellipse.a ** 2 * tline.v.y ** 2 -
              tline.v.y ** 2 * tline.p1.x ** 2 +
              2 * tline.v.x * tline.v.y * tline.p1.x * tline.p1.y -
              tline.v.x ** 2 * tline.p1.y ** 2
          ) *
            ellipse.a *
            ellipse.b
        ) /
        (ellipse.b ** 2 * tline.v.x ** 2 + ellipse.a ** 2 * tline.v.y ** 2);
      let t2 =
        -(
          ellipse.b ** 2 * tline.v.x * tline.p1.x +
          ellipse.a ** 2 * tline.v.y * tline.p1.y -
          Math.sqrt(
            ellipse.b ** 2 * tline.v.x ** 2 +
              ellipse.a ** 2 * tline.v.y ** 2 -
              tline.v.y ** 2 * tline.p1.x ** 2 +
              2 * tline.v.x * tline.v.y * tline.p1.x * tline.p1.y -
              tline.v.x ** 2 * tline.p1.y ** 2
          ) *
            ellipse.a *
            ellipse.b
        ) /
        (ellipse.b ** 2 * tline.v.x ** 2 + ellipse.a ** 2 * tline.v.y ** 2);
      if ((t1 >= 0 && t1 <= tline.t) || (t2 >= 0 && t2 <= tline.t)) {
        return true;
      }
    }

    // If the line doesn't intersect the ellipse, then we just check its endpoints - it could still be entirely within the ellipse
    // Technically this could equivalently be && instead of || but no reason not to play it safe
    return (
      Hitbox.collideCoordEllipse({ type: "coords", p: line.p1 }, ellipse) ||
      Hitbox.collideCoordEllipse({ type: "coords", p: line.p2 }, ellipse)
    );
  }

  static collidePolyPoly(poly1, poly2) {
    for (let i = 0; i < poly1.ps.length; i++) {
      let line1 = Hitbox.makeLine(
        poly1.ps[i],
        poly1.ps[(i + 1) % poly1.ps.length]
      );
      for (let j = 0; j < poly2.ps.length; j++) {
        let line2 = Hitbox.makeLine(
          poly2.ps[j],
          poly2.ps[(j + 1) % poly2.ps.length]
        );
        if (Hitbox.collideLineLine(line1, line2)) {
          return true;
        }
      }
    }

    return (
      [...Array(poly1.ps.length).keys()].reduce(
        (l, i) =>
          l || Hitbox.collideCoordPoly({ type: "coord", p: poly1.ps[i] }, poly2)
      ) ||
      [...Array(poly2.ps.length).keys()].reduce(
        (l, i) =>
          l || Hitbox.collideCoordPoly({ type: "coord", p: poly2.ps[i] }, poly1)
      )
    );
  }

  static collidePolyEllipse(poly, ellipse) {
    for (let i = 0; i < poly.ps.length; i++) {
      let line = Hitbox.makeLine(poly.ps[i], poly.ps[(i + 1) % poly.ps.length]);

      if (Hitbox.collideLineEllipse(line, ellipse)) {
        return true;
      }
    }

    return (
      [...Array(poly.ps.length).keys()].reduce(
        (l, i) =>
          l ||
          Hitbox.collideCoordEllipse({ type: "coord", p: poly.ps[i] }, ellipse)
      ) || Hitbox.collideCoordPoly({ type: "coord", p: ellipse.center }, poly)
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
        .transformPoint(new DOMPoint(ellipse.center.x, ellipse.center.y));
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
      Hitbox.collideCoordEllipse(
        { type: "coord", p: ellipse1.center },
        ellipse2
      ) ||
      Hitbox.collideCoordEllipse(
        { type: "coord", p: ellipse2.center },
        ellipse1
      )
    ) {
      return true;
    }

    // Convert both ellipses to bivariate form, then check if they intersect
    let elps1 = bivariateForm(ellipse1);
    let elps2 = bivariateForm(ellipse2);
    return conicsIntersect(elps1, elps2);
  }
}
