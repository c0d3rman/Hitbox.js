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
    for (var i = 0; i < this.components.length; i++) {
      for (var j = 0; j < other.components.length; j++) {
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

  line(x1, y1, x2, y2) {
    if (this.draw) {
      line(x1, y1, x2, y2);
    }
    [x1, y1] = Hitbox.transformPoint(x1, y1);
    [x2, y2] = Hitbox.transformPoint(x2, y2);
    this.components.push(Hitbox.makeLine(x1, y1, x2, y2));
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

  static collideEllipseEllipse(ellipse1, ellipse2) {
    // First check if one ellipse is within the other
    if (
      Hitbox.collideCoordEllipse({ x: ellipse1.x, y: ellipse1.y }, ellipse2) ||
      Hitbox.collideCoordEllipse({ x: ellipse2.x, y: ellipse2.y }, ellipse1)
    ) {
      return true;
    }

    /*
    1. Move plane so ellipse1 is centered at the origin
    2. Rotate plane so ellipse1 is flat (major axis parallel to x axis)
    3. Transform plane so ellipse1 is a circle
    (order of operations is reversed due to how matrix math works)
    */
    let mat1 = new DOMMatrix()
      .scale(ellipse1.b / ellipse1.a, 1)
      .rotate((180 / PI) * -ellipse1.angle)
      .translate(-ellipse1.x, -ellipse1.y);
    ellipse1 = Hitbox.transformEllipse(ellipse1, mat1);
    ellipse2 = Hitbox.transformEllipse(ellipse2, mat1);

    /*
    4. Move plane so ellipse2 is centered at the origin
    5. Rotate plane so ellipse2 is flat (major axis parallel to x axis)
    */
    let mat2 = new DOMMatrix()
      .rotate((180 / PI) * -ellipse2.angle)
      .translate(-ellipse2.x, -ellipse2.y);
    ellipse1 = Hitbox.transformEllipse(ellipse1, mat2);
    ellipse2 = Hitbox.transformEllipse(ellipse2, mat2);

    // With the parametric form of ellipse2 as (a cos t, b sin t),
    // and the center of circle ellipse1 as (x, y),
    // Use Newton approximation to find the t that minimizes
    // (a * cos(t) - x)^2 + (b * sin(t) - y)^2
    // This is the closest point on ellipse2 to the center of the circle.
    // I really wanted a closed form of this but turns out that even though it exists, it's a very big and complex equation.

    let f = (t) =>
      2 * (ellipse2.b * sin(t) - ellipse1.y) * ellipse2.b * cos(t) -
      2 * (ellipse2.a * cos(t) - ellipse1.x) * ellipse2.a * sin(t);

    let g = (t) =>
      2 * ellipse2.b ** 2 * cos(t) ** 2 +
      2 * ellipse2.a ** 2 * sin(t) ** 2 -
      2 * (ellipse2.a * cos(t) - ellipse1.x) * ellipse2.a * cos(t) -
      2 * (ellipse2.b * sin(t) - ellipse1.y) * ellipse2.b * sin(t);

    let t0 = Math.atan2(ellipse1.y, ellipse1.x);
    let min_error = 0.1;
    let max_iter = 100;
    let step = 0;

    while (abs(f(t0)) > min_error) {
      step++;

      if (g(t0) == 0) {
        console.log("Mathematical Error");
        return false;
      }

      t0 = t0 - f(t0) / g(t0);

      if (step > max_iter) {
        console.log("Not Convergent");
        return false;
      }
    }

    // Now that we know the closest point on ellipse2 to the center of the circle ellipse1, it's just a question of distance
    return (
      Hitbox.distance(
        ellipse2.a * cos(t0),
        ellipse2.b * sin(t0),
        ellipse1.x,
        ellipse1.y
      ) <= ellipse1.a
    );
  }
}