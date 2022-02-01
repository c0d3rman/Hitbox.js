let EPSILON = 0.00001;

// Create a canvas element for internal use in detecting collisions
_hitbox_canvas = document.createElement("canvas");

// Mode constants
DRAW_ON = "draw_on";
DRAW_OFF = "draw_off";
FILL_ON = "fill_on";
FILL_OFF = "fill_off";
RESPECT_STROKE_WEIGHT = "respect_stroke_weight";
IGNORE_STROKE_WEIGHT = "ignore_stroke_weight";

// For heuristic checking of whether components are close enough to each other to be worth checking for a collision.
// Saves a lot of time performing expensive ellipse computations.
class HitboxBoundingCircle {
  constructor(center, radius) {
    this.center = center;
    this.radius = radius;
  }

  collidesWith(other) {
    return (
      this.center.dist(other.center) <= this.radius + other.radius + EPSILON
    );
  }

  contains(p) {
    return p.dist(this.center) <= this.radius + EPSILON;
  }

  // Returns a shuffled copy of an array.
  // Adapted from https://stackoverflow.com/a/12646864/2674563
  static shuffle(array) {
    array = [...array];
    for (let i = array.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [array[i], array[j]] = [array[j], array[i]];
    }
    return array;
  }

  // Returns the smallest circle that encloses all the given points.
  // Uses Welzl's algorithm.
  // If you know any points that  must be on the circle's edge, you can pass them in as R.
  static fromPoints(P, R = [], didShuffle = false) {
    // Only shuffle P once
    if (!didShuffle) {
      P = [...P]; // Make a copy of P to not interfere with the caller
      HitboxBoundingCircle.shuffle(P);
    }

    // Base case
    if (P.length == 0 || R.length == 3) {
      if (R.length == 0) {
        return new HitboxBoundingCircle(createVector(0, 0), -1); // Dummy circle with negative radius that contains nothing
      }

      if (R.length == 1) {
        return new HitboxBoundingCircle(R[0], 0);
      }

      if (R.length == 2) {
        return new HitboxBoundingCircle(
          p5.Vector.add(R[0], R[1]).div(2),
          R[0].dist(R[1]) / 2
        );
      }

      // 3-point case
      // Math from wikipedia.org/wiki/Circumscribed_circle#Circumcenter_coordinates
      let o = createVector(
        (Math.min(R.map((p) => p.x)) + Math.max(R.map((p) => p.x))) / 2,
        (Math.min(R.map((p) => p.y)) + Math.max(R.map((p) => p.y))) / 2
      );
      let a = p5.Vector.sub(R[0], o);
      let b = p5.Vector.sub(R[1], o);
      let c = p5.Vector.sub(R[2], o);
      let d = (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) * 2;

      if (d == 0) {
        // The three points form a line, so get the two outer ones and recurse to the two-point case
        return HitboxBoundingCircle.fromPoints(
          [],
          [
            R.reduce((p1, p2) => (p1.x < p2.x ? p1 : p2)),
            R.reduce((p1, p2) => (p1.x > p2.x ? p1 : p2)),
          ],
          true
        );
      }

      o.add(
        createVector(
          (a.magSq() * (b.y - c.y) +
            b.magSq() * (c.y - a.y) +
            c.magSq() * (a.y - b.y)) /
            d,
          (a.magSq() * (c.x - b.x) +
            b.magSq() * (a.x - c.x) +
            c.magSq() * (b.x - a.x)) /
            d
        )
      );

      return new HitboxBoundingCircle(o, o.dist(a));
    }

    // Get the first point (which is random because we shuffled P)
    let p = P[0];
    // Recursively find the smallest circle which includes all of the other points
    let D = HitboxBoundingCircle.fromPoints(P.slice(1), R, true);
    // If p is in said circle already, then that's our cirlce
    if (D.contains(p)) {
      return D;
    }

    // If not, p must be on the boundary of the true circle - recurse with that assumption in place
    return HitboxBoundingCircle.fromPoints(P.slice(1), [...R, p], true);
  }
}

class HitboxCoord {
  constructor(p) {
    this.type = "Coord";
    this.includeFill = false; // No fill
    this.p = p;

    // For checking if the we're entirely inside a filled shape
    this.representativePoint = p;
  }

  draw() {
    push();
    strokeWeight(5);
    point(this.p.x, this.p.y);
    pop();
  }

  smear(v) {
    return [new HitboxLine(this.p, p5.Vector.add(this.p, v))];
  }
}

class HitboxLine {
  constructor(p1, p2) {
    this.type = "Line";
    this.includeFill = false; // No fill

    // Endpoints of the line
    this.p1 = p1;
    this.p2 = p2;

    // Parametric representation of the line as p = p1 + vt
    this.v = p5.Vector.sub(p2, p1).normalize();
    this.t = p2.dist(p1);

    // For checking if the we're entirely inside a filled shape
    this.representativePoint = p1;

    this.boundingCircle = HitboxBoundingCircle.fromPoints([], [p1, p2], true); // P is "shuffled" because it's empty
  }

  // Check which side of the line a point is on
  sideOf(p) {
    if (Hitbox.collideCoordLine(new HitboxCoord(p), this)) {
      return 0;
    }

    // We draw a vector between the line's starting point and p, and then cross it with the line's direction vector
    return Math.sign(p5.Vector.sub(this.p1, p).cross(this.v).z);
  }

  draw() {
    line(this.p1.x, this.p1.y, this.p2.x, this.p2.y);
  }

  smear(v) {
    let cosAngle = v.dot(this.v) / (v.mag() * this.v.mag());
    // Parallel case - θ = 0, cos(θ) = 1
    if (Math.abs(cosAngle) < EPSILON) {
      return this.p1.smear(v);
    }
    // Antiparallel case - θ = π, cos(θ) = -1
    else if (Math.abs(cosAngle + 1) < EPSILON) {
      return this.p2.smear(v);
    }

    // Normal case
    return [
      new HitboxPoly(
        [
          this.p1,
          this.p2,
          p5.Vector.add(v, this.p2),
          p5.Vector.add(v, this.p1),
        ],
        true
      ),
    ];
  }
}

class HitboxPoly {
  constructor(ps, includeFill) {
    this.type = "Poly";

    // Vertices of the polygon
    this.ps = ps;

    // Should this be an outline-only polygon, or should it include the inside
    this.includeFill = includeFill;
    this.representativePoint = ps[0];

    // Make the lines
    this.lines = [];
    for (let i = 0; i < ps.length; i++) {
      this.lines.push(new HitboxLine(ps[i], ps[(i + 1) % ps.length]));
    }

    // Make a canvas path for use in coordInsidePoly
    this.path = new Path2D();
    this.path.moveTo(this.ps[0].x, this.ps[0].y);
    for (let i = 1; i < this.ps.length; i++) {
      this.path.lineTo(ps[i].x, ps[i].y);
    }
    this.path.closePath();

    // Make a bounding circle
    this.boundingCircle = HitboxBoundingCircle.fromPoints(ps);
  }

  draw() {
    beginShape();
    for (const p of this.ps) {
      vertex(p.x, p.y);
    }
    endShape(CLOSE);
  }

  smear(v) {
    // Find the two outermost points perpendicular to the direction of the smear
    let perpV = createVector(v.y, -v.x).normalize();
    let projectedYs = this.ps.map((p) => p.dot(perpV));
    let i_min = this.ps.reduce(
      (i_prev, _, i) => (projectedYs[i] < projectedYs[i_prev] ? i : i_prev),
      0
    );
    let i_max = this.ps.reduce(
      (i_prev, _, i) => (projectedYs[i] > projectedYs[i_prev] ? i : i_prev),
      0
    );

    // Draw two smear-polygons with the polygon split into two halves using said points
    let poly = Hitbox.transform(this, new DOMMatrix().translate(v.x, v.y));
    let ps1 = [];
    let ps2 = [];
    for (let i = i_min; i != i_max; i = (i + 1) % this.ps.length) {
      ps1.push(this.ps[i]);
      ps2.push(poly.ps[i]);
    }
    ps1.push(this.ps[i_max]);
    ps2.push(poly.ps[i_max]);
    for (let i = i_max; i != i_min; i = (i + 1) % this.ps.length) {
      ps1.push(poly.ps[i]);
      ps2.push(this.ps[i]);
    }
    ps1.push(poly.ps[i_min]);
    ps2.push(this.ps[i_min]);

    return [new HitboxPoly(ps1, true), new HitboxPoly(ps2, true)];
  }
}

class HitboxEllipse {
  // The ellipse is defined by two conjugate diameters.
  // This is because conjugate diameters are resilient to translation, rotation, shearing, and the like.
  // We specify the cojugate diameters by giving one endpoint of each, as well as the center of the ellipse.
  // This math is from wikipedia: https://en.wikipedia.org/wiki/Ellipse#General_ellipse_2
  constructor(center, dia1, dia2, includeFill) {
    this.type = "Ellipse";

    this.center = center;
    this.dia1 = dia1;
    this.dia2 = dia2;

    this.includeFill = includeFill;
    this.representativePoint = center;

    let vf1 = p5.Vector.sub(dia1, center);
    let vf2 = p5.Vector.sub(dia2, center);

    // Get the t of the first conjugate diameter
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
    this.a = center.dist(v1);
    this.b = center.dist(v2);

    // Make sure a is the major and b is the minor
    if (this.b > this.a) {
      [this.a, this.b] = [this.b, this.a];
      [v1, v2] = [v2, v1];
    }

    // Calculate angle of rotation
    this.angle = Math.atan2(v1.y - center.y, v1.x - center.x);

    // Calculate foci using semiaxes
    let c = Math.sqrt(this.a ** 2 - this.b ** 2);
    // Thankfully fromAngle is unaffected by angleMode
    this.f1 = p5.Vector.add(center, p5.Vector.fromAngle(this.angle, c));
    this.f2 = p5.Vector.sub(center, p5.Vector.fromAngle(this.angle, c));

    // r is the sum of distances between any point on the ellipse and the two foci
    // It's also just the major axis length
    this.r = this.a * 2;

    this.boundingCircle = new HitboxBoundingCircle(center, this.a);
  }

  draw() {
    push();
    translate(this.center.x, this.center.y);
    rotate(this.angle);
    ellipse(0, 0, this.a * 2, this.b * 2);
    pop();
  }
}

class HitboxArc {
  // The arc is defined by an ellipse and two endpoints on it.
  // This is because points are resilient to translation, rotation, shearing, and the like (unlike angles).
  constructor(ellipse, startV, stopV, mode) {
    this.type = "Arc";

    this.ellipse = ellipse;
    this.startV = startV;
    this.stopV = stopV;
    this.mode = mode;

    this.includeFill = ellipse.includeFill;
    this.representativePoint = startV;

    this.start = Math.atan2(
      startV.y - ellipse.center.y,
      startV.x - ellipse.center.x
    );
    this.stop = Math.atan2(
      stopV.y - ellipse.center.y,
      stopV.x - ellipse.center.x
    );

    // Deal with the non-curved edges of the arc, based on the arc mode
    this.lines = [];
    if (mode == p5.prototype.PIE || (this.includeFill && mode == null)) {
      this.lines.push(new HitboxLine(this.ellipse.center, this.startV));
      this.lines.push(new HitboxLine(this.ellipse.center, this.stopV));
    } else if (
      mode == p5.prototype.CHORD ||
      (this.includeFill && mode == p5.prototype.OPEN)
    ) {
      this.lines.push(new HitboxLine(this.startV, this.stopV));
    }

    // Not the optimal bounding circle for small arcs, but who cares
    this.boundingCircle = new HitboxBoundingCircle(ellipse.center, ellipse.a);
  }

  draw() {
    push();
    translate(this.ellipse.center.x, this.ellipse.center.y);
    rotate(this.ellipse.angle);
    arc(
      0,
      0,
      this.ellipse.a * 2,
      this.ellipse.b * 2,
      this.start - this.ellipse.angle,
      this.stop - this.ellipse.angle,
      this.mode
    );
    pop();
  }
}

// static makeLine() {

//     if (!this.includeStroke || !makeStrokePolygon) {
//       return l;
//     } else {
//       let perpV = createVector(l.v.y, -l.v.x) // perpendicular vector to line
//         .normalize()
//         .mult(drawingContext.lineWidth / 2);

//       return {
//         type: "poly",
//         ps: [
//           p5.Vector.add(l.p1, perpV),
//           p5.Vector.sub(l.p1, perpV),
//           p5.Vector.sub(l.p2, perpV),
//           p5.Vector.add(l.p2, perpV),
//         ],
//       };
//     }
//   }

class Hitbox {
  constructor(...args) {
    // Default values
    this.drawMode = DRAW_OFF;
    this.fillMode = FILL_ON;
    this.strokeWeightMode = IGNORE_STROKE_WEIGHT;

    // Update values if any given
    let { drawMode, fillMode, strokeWeightMode } = this.handleModeArgs(args);
    this.drawMode = drawMode;
    this.fillMode = fillMode;
    this.strokeWeightMode = strokeWeightMode;

    if (args.length > 0) {
      throw "Illegal arguments passed to Hitbox constructor";
    }

    this.components = [];

    this.collisionFnMap = {
      coordInsidePoly: Hitbox.coordInsidePoly,
      coordInsideEllipse: Hitbox.coordInsideEllipse,
      coordInsideArc: Hitbox.coordInsideArc,
      collideCoordCoord: Hitbox.collideCoordCoord,
      collideCoordLine: Hitbox.collideCoordLine,
      collideCoordPoly: Hitbox.collideCoordPoly,
      collideCoordEllipse: Hitbox.collideCoordEllipse,
      collideCoordArc: Hitbox.collideCoordArc,
      collideLineLine: Hitbox.collideLineLine,
      collideLinePoly: Hitbox.collideLinePoly,
      collideLineEllipse: Hitbox.collideLineEllipse,
      collideLineArc: Hitbox.collideLineArc,
      collidePolyPoly: Hitbox.collidePolyPoly,
      collidePolyEllipse: Hitbox.collidePolyEllipse,
      collidePolyArc: Hitbox.collidePolyArc,
      collideEllipseEllipse: Hitbox.collideEllipseEllipse,
      collideEllipseArc: Hitbox.collideEllipseArc,
      collideArcArc: Hitbox.collideArcArc,
    };
  }

  // Public functions - collisions & other
  // =====================================

  collidesWith(other) {
    for (let A of this.components) {
      for (let B of other.components) {
        // First, check bounding circles if they exist. If the bounding circles don't touch, no reason to keep going.
        if (
          "boundingCircle" in A &&
          "boundingCircle" in B &&
          !A.boundingCircle.collidesWith(B.boundingCircle)
        ) {
          continue;
        }

        // Now, handle fill: if component A includes its fill, check if component B is entirely within component A (and vice versa) by testing a representative point
        if (
          A.includeFill &&
          this.collisionFnMap["coordInside" + A.type](
            new HitboxCoord(B.representativePoint),
            A
          )
        ) {
          return true;
        }
        if (
          B.includeFill &&
          this.collisionFnMap["coordInside" + B.type](
            new HitboxCoord(A.representativePoint),
            B
          )
        ) {
          return true;
        }

        // Now check the actual collision function, swapping A and B if necessary
        let key = "collide" + A.type + B.type;
        if (key in this.collisionFnMap) {
          if (this.collisionFnMap[key](A, B)) {
            return true;
          }
        } else {
          key = "collide" + B.type + A.type;
          if (this.collisionFnMap[key](B, A)) {
            return true;
          }
        }
      }
    }
    return false;
  }

  collidesWithPoint(x, y) {
    return this.collidesWith({
      components: [createVector(x, y)],
    });
  }

  // Draws the hitbox components as they actually are.
  // For debugging purposes.
  drawDebug() {
    push();
    resetMatrix();
    for (const c of this.components) {
      if (c.includeFill) {
        fill(color(255, 0, 0, 100));
        noStroke();
      } else {
        noFill();
        stroke(255, 0, 0);
      }
      c.draw();
    }
    pop();
    return this;
  }

  // Turn drawing on and off.
  // Valid values:
  // DRAW_OFF - don't draw (default)
  // DRAW_ON - draw like p5 does in addition to creating the hitbox
  setDrawMode(drawMode) {
    if (![DRAW_ON, DRAW_OFF].includes(drawMode)) {
      throw "Invalid draw mode";
    }
    this.drawMode = drawMode;
    return this;
  }

  // Change whether to use shapes as just their outlines,
  // or whether to include their inside area as part of the hitbox.
  // Valid values:
  // FILL_ON - include (default)
  // FILL_OFF - don't include
  setFillMode(fillMode) {
    if (![FILL_ON, FILL_OFF].includes(fillMode)) {
      throw "Invalid fill mode";
    }
    this.includeStroke = includeStroke;
    return this;
  }

  // Change whether to treat the stroke of a shape as a one-pixel boundary,
  // or whether to take its stroke weight into account (e.g. turn a point into a circle).
  // Valid values:
  // IGNORE_STROKE_WEIGHT - treat as one-pixel boundary (default)
  // RESPECT_STROKE_WEIGHT - take stroke weight into account
  setStrokeWeightMode(strokeWeightMode) {
    if (
      ![RESPECT_STROKE_WEIGHT, IGNORE_STROKE_WEIGHT].includes(strokeWeightMode)
    ) {
      throw "Invalid stroke weight mode";
    }
    this.strokeWeightMode = strokeWeightMode;
    return this;
  }

  // Public functions - component creation
  // =====================================
  
  startCapture() {
    if (Hitbox.prototype.currentCapturingHitbox != null) {
      throw "Can't start capture because there is already a capture in progress"
    }
    if (this.drawMode == DRAW_ON) {
      console.log("Hitbox warning: draw mode will be set to DRAW_OFF during capturing to avoid infinite looping.")
    }
    Hitbox.prototype.currentCapturingHitbox = this
  }
    
  stopCapture() {
    if (Hitbox.prototype.currentCapturingHitbox !== this) {
      throw "This hitbox is not in the process of a capture"
    }
    Hitbox.prototype.currentCapturingHitbox = null
  }

  point(...args) {
    let { drawMode, strokeWeightMode } = this.handleModeArgs(args); // Ignore fillMode

    if (drawMode == DRAW_ON) {
      point(...args);
    }

    let p;
    if (args.length === 1 && args[0] instanceof p5.Vector) {
      p = args[0];
    } else {
      p = createVector(args[0], args[1]);
    }

    if (strokeWeightMode == IGNORE_STROKE_WEIGHT) {
      this.components.push(Hitbox.transform(new HitboxCoord(p)));
    } else {
      this.components.push(
        Hitbox.transform(
          new HitboxEllipse(
            p,
            createVector(p.x + drawingContext.lineWidth / 2, p.y),
            createVector(p.x, p.y + drawingContext.lineWidth / 2),
            true
          )
        )
      );
    }

    return this;
  }

  line(...args) {
    let { drawMode, strokeWeightMode } = this.handleModeArgs(args); // Ignore fillMode

    let [x1, y1, x2, y2] = args;

    if (drawMode == DRAW_ON) {
      line(x1, y1, x2, y2);
    }

    let l = Hitbox.transform(
      new HitboxLine(createVector(x1, y1), createVector(x2, y2))
    );

    if (strokeWeightMode == IGNORE_STROKE_WEIGHT) {
      this.components.push(l);
    } else {
      let perpV = createVector(l.v.y, -l.v.x) // perpendicular vector to line
        .normalize()
        .mult(drawingContext.lineWidth / 2);
      this.components.push(
        Hitbox.transform(
          new HitboxPoly([
            p5.Vector.add(l.p1, perpV),
            p5.Vector.sub(l.p1, perpV),
            p5.Vector.sub(l.p2, perpV),
            p5.Vector.add(l.p2, perpV),
          ])
        )
      );

      if (includeStrokeCap) {
        switch (drawingContext.lineCap) {
          case p5.prototype.ROUND:

          case p5.prototype.SQUARE:
            // Do nothing
            break;
          case p5.prototype.PROJECT:
            break;
        }
      }
    }

    return this;
  }

  poly(...args) {
    let { drawMode, fillMode, strokeWeightMode } = this.handleModeArgs(args);

    let ps = args[0];

    if (drawMode == DRAW_ON) {
      beginShape();
      for (const p of ps) {
        vertex(p.x, p.y);
      }
      endShape(CLOSE);
    }

    this.components.push(
      Hitbox.transform(new HitboxPoly(ps, fillMode == FILL_ON))
    );
  }

  rect(...args) {
    let { drawMode, fillMode, strokeWeightMode } = this.handleModeArgs(args);

    if (drawMode == DRAW_ON) {
      rect(...args);
    }

    // Add in the height if none provided
    if (args.length == 3) {
      args.splice(3, 0, args[2]);
    }

    let { x, y, w, h } = Hitbox.modeAdjust(
      args[0],
      args[1],
      args[2],
      args[3],
      window._renderer._rectMode
    );

    // p5 allows negative heights/widths
    w = Math.abs(w);
    h = Math.abs(h);

    // Rounded corner rect
    if (args.length > 4) {
      // Fill all 4 rounded-corner params if some are not provided
      while (args.length < 8) {
        args.push(args[args.length - 1]);
      }

      // Clip the radii to at most half the width or the height
      let [tl, tr, br, bl] = args
        .slice(4)
        .map((r) => Math.min(r, w / 2, h / 2));

      // Draw an octagon for the non-curvy part of the square (chopped-off corners)
      let ps = [
        createVector(x, y + tl),
        createVector(x + tl, y),
        createVector(x + w - tr, y),
        createVector(x + w, y + tr),
        createVector(x + w, y + h - br),
        createVector(x + w - br, y + h),
        createVector(x + bl, y + h),
        createVector(x, y + h - bl),
      ];
      this.poly(ps, DRAW_OFF, fillMode, strokeWeightMode);

      // Draw the four rounded corner arcs
      this.arc(
        x + tl,
        y + tl,
        tl * 2,
        tl * 2,
        PI,
        PI * (3 / 2),
        DRAW_OFF,
        fillMode,
        strokeWeightMode
      );
      this.arc(
        x + w - tr,
        y + tr,
        tr * 2,
        tr * 2,
        PI * (3 / 2),
        0,
        DRAW_OFF,
        fillMode,
        strokeWeightMode
      );
      this.arc(
        x + w - br,
        y + h - br,
        br * 2,
        br * 2,
        0,
        PI / 2,
        DRAW_OFF,
        fillMode,
        strokeWeightMode
      );
      this.arc(
        x + bl,
        y + h - bl,
        bl * 2,
        bl * 2,
        PI / 2,
        PI,
        PIE,
        DRAW_OFF,
        fillMode,
        strokeWeightMode
      );

      return this;
    }

    // Regular rect
    let ps = [
      createVector(x, y),
      createVector(x + w, y),
      createVector(x + w, y + h),
      createVector(x, y + h),
    ];
    return this.poly(ps, DRAW_OFF, fillMode, strokeWeightMode);
  }

  triangle(...args) {
    let { drawMode, fillMode, strokeWeightMode } = this.handleModeArgs(args);

    let [x1, y1, x2, y2, x3, y3] = args;

    if (drawMode == DRAW_ON) {
      triangle(x1, y1, x2, y2, x3, y3);
    }

    let ps = [createVector(x1, y1), createVector(x2, y2), createVector(x3, y3)];

    return this.poly(ps, DRAW_OFF, fillMode, strokeWeightMode);
  }

  quad(...args) {
    let { drawMode, fillMode, strokeWeightMode } = this.handleModeArgs(args);

    let [x1, y1, x2, y2, x3, y3, x4, y4] = args;

    if (drawMode == DRAW_ON) {
      quad(x1, y1, x2, y2, x3, y3, x4, y4);
    }

    let ps = [
      createVector(x1, y1),
      createVector(x2, y2),
      createVector(x3, y3),
      createVector(x4, y4),
    ];

    return this.poly(ps, DRAW_OFF, fillMode, strokeWeightMode);
  }

  ellipse(...args) {
    let { drawMode, fillMode, strokeWeightMode } = this.handleModeArgs(args);

    if (drawMode == DRAW_ON) {
      ellipse(...args);
    }

    // Support for no-height syntax
    if (args.length == 3) {
      args.push(args[2]);
    }

    let { x, y, w, h } = Hitbox.modeAdjust(
      args[0],
      args[1],
      args[2],
      args[3],
      window._renderer._ellipseMode
    );

    // we want things in CENTER mode, not CORNER
    x += w / 2;
    y += h / 2;

    // Handle stroke
    // Technically the resulting shape should be an octic, but p5's renderer just draws an ellipse so whatever
    if (this.includeStroke) {
      w += drawingContext.lineWidth;
      h += drawingContext.lineWidth;
    }

    // Get conjugate diameters by transforming the original vertices of the ellipse
    this.components.push(
      Hitbox.transform(
        new HitboxEllipse(
          createVector(x, y),
          createVector(x + w / 2, y),
          createVector(x, y + h / 2),
          fillMode == FILL_ON
        )
      )
    );

    return this;
  }

  arc(...args) {
    let { drawMode, fillMode, strokeWeightMode } = this.handleModeArgs(args);

    let [x, y, w, h, start, stop] = args.slice(0, 6);
    let mode = args.length > 6 ? args[6] : null; // If no mode arg, mode is null

    if (drawMode == DRAW_ON) {
      arc(...args);
    }

    // Normalize angles just like p5 does
    const angles = p5.prototype._normalizeArcAngles(start, stop, w, h, true);

    if (angles.correspondToSamePoint) {
      // Just use an ellipse hitbox in this case
      return this.ellipse(x, y, w, h, DRAW_OFF, fillMode, strokeWeightMode);
    }

    ({ x, y, w, h } = Hitbox.modeAdjust(
      x,
      y,
      w,
      h,
      window._renderer._ellipseMode
    ));

    // we want things in CENTER mode, not CORNER
    x += w / 2;
    y += h / 2;

    // Handle stroke
    // Technically the resulting shape should be an octic, but p5's renderer just draws an ellipse so whatever
    if (strokeWeightMode == RESPECT_STROKE_WEIGHT) {
      w += drawingContext.lineWidth;
      h += drawingContext.lineWidth;
    }

    // Get conjugate diameters by transforming the original vertices of the ellipse
    this.components.push(
      Hitbox.transform(
        new HitboxArc(
          new HitboxEllipse(
            createVector(x, y),
            createVector(x + w / 2, y),
            createVector(x, y + h / 2),
            fillMode == FILL_ON
          ),
          createVector(
            x + (w / 2) * cos(angles.start),
            y + (h / 2) * sin(angles.start)
          ),
          createVector(
            x + (w / 2) * cos(angles.stop),
            y + (h / 2) * sin(angles.stop)
          ),
          mode
        )
      )
    );

    return this;
  }

  square(...args) {
    // Add in the height
    args.splice(3, 0, args[2]);

    return this.rect(...args);
  }

  circle(...args) {
    // Add in the height
    args.splice(3, 0, args[2]);

    return this.ellipse(...args);
  }

  // Private helpers
  // ===============

  handleModeArgs(args) {
    let out = {};

    for (let i = args.length - 1; i >= 0; i--) {
      let arg = args[i];
      if ([DRAW_ON, DRAW_OFF].includes(arg)) {
        if ("drawMode" in out) {
          throw "Can't use more than one draw mode argument";
        }
        out.drawMode = arg;
      } else if ([FILL_ON, FILL_OFF].includes(arg)) {
        if ("fillMode" in out) {
          throw "Can't use more than one fill mode argument";
        }
        out.fillMode = arg;
      } else if ([RESPECT_STROKE_WEIGHT, IGNORE_STROKE_WEIGHT].includes(arg)) {
        if ("strokeWeightMode" in out) {
          throw "Can't use more than one stroke weight mode argument";
        }
        out.strokeWeightMode = arg;
      } else {
        continue;
      }
      args.splice(i, 1);
    }

    // If no value was given, use current default
    if (!("drawMode" in out)) {
      out.drawMode = this.drawMode;
    }
    if (!("fillMode" in out)) {
      out.fillMode = this.fillMode;
    }
    if (!("strokeWeightMode" in out)) {
      out.strokeWeightMode = this.strokeWeightMode;
    }

    return out;
  }

  // Duplicated function from inside p5. Handles rectMode and ellipseMode.
  // Modified to also handle allowing negative widths and heights.
  static modeAdjust(a, b, c, d, mode) {
    let out;
    if (mode === p5.prototype.CORNER) {
      out = { x: a, y: b, w: c, h: d };
    } else if (mode === p5.prototype.CORNERS) {
      out = { x: a, y: b, w: c - a, h: d - b };
    } else if (mode === p5.prototype.RADIUS) {
      out = { x: a - c, y: b - d, w: 2 * c, h: 2 * d };
    } else if (mode === p5.prototype.CENTER) {
      out = { x: a - c * 0.5, y: b - d * 0.5, w: c, h: d };
    } else {
      throw "Unknown mode";
    }

    if (out.w < 0) {
      out.w = Math.abs(out.w);
      out.x -= out.w;
    }
    if (out.h < 0) {
      out.h = Math.abs(out.h);
      out.y -= out.h;
    }

    return out;
  }

  // Check if angle is an angle between the start and stop angles,
  // handling periodicity.
  // All angles must be in radians.
  static isAngleBetween(angle, start, stop) {
    // Make sure all angles are in [0, 2π)
    angle -= TWO_PI * Math.floor(angle / TWO_PI);
    start -= TWO_PI * Math.floor(start / TWO_PI);
    stop -= TWO_PI * Math.floor(stop / TWO_PI);

    if (start < stop) {
      return angle >= start && angle <= stop;
    } else {
      return angle >= start || angle <= stop;
    }
  }

  // Given a cubic of the form y³ + b*y² + c*y + d = 0,
  // find one real root. (There is always at least one.)
  // This follows the formula from https://math.vanderbilt.edu/schectex/courses/cubic/
  // Used for ellipse-ellipse collisions.
  static solveCubic(b, c, d) {
    // We need to do the math in a way that can handle complex numbers.
    // Even though the coefficients and roots are both real, the process
    // of solving can involve intermediate steps with complex numbers.
    // So we use mathjs - this line pulls out relevant functions from it
    // so we can use shorthand versions of them.
    let {
      add,
      subtract: sub,
      multiply: mul,
      divide: div,
      pow,
      sqrt,
      cbrt,
    } = math;

    // Calculate the real root.
    // (The other two roots of the cubic might also be real, but we don't care.)
    let p = add(div(pow(b, 3), -27), mul(b, c, 1 / 6), div(d, -2));
    let q = sqrt(add(pow(p, 2), pow(sub(div(c, 3), div(pow(b, 2), 9)), 3)));
    let y_1 = add(cbrt(add(p, q)), cbrt(sub(p, q)), div(b, -3));

    // Make sure the result is returned as a real number instead of a complex
    // number with a zero imaginary component.
    if (!math.isNumeric(y_1)) {
      if (y_1.im != 0) {
        throw "Something went very wrong finding the real root of a cubic.";
      }
      y_1 = y_1.re;
    }

    return y_1;
  }

  // Given a quartic of the form z_4*y⁴ + z_3*y³ + z_2*y² + z_1*y + z_0 = 0,
  // find all of its roots.
  // This follows the solution from https://mathworld.wolfram.com/QuarticEquation.html - particularly equations 34 and 36-42.
  // Used for ellipse-ellipse collisions.
  static solveQuartic(z_4, z_3, z_2, z_1, z_0) {
    // Normalize coefficients by z_4 to make z_4 = 1
    let [a_3, a_2, a_1, a_0] = [z_3 / z_4, z_2 / z_4, z_1 / z_4, z_0 / z_4];

    // Get the resolvent cubic's coefficients (equation 34 from Wolfram)
    // a always equals 1
    let b = -a_2;
    let c = a_1 * a_3 - 4 * a_0;
    let d = 4 * a_2 * a_0 - a_1 ** 2 - a_3 ** 2 * a_0;

    // Find a real root of the cubic
    let y_1 = Hitbox.solveCubic(b, c, d);

    // Use equations 36-42 from Wolfram to find the four roots
    let R = sqrt((1 / 4) * a_3 ** 2 - a_2 + y_1);
    let D, E;
    if (R != 0) {
      D = sqrt(
        (3 / 4) * a_3 ** 2 -
          R ** 2 -
          2 * a_2 +
          (4 * a_3 * a_2 - 8 * a_1 - a_3 ** 3) / (4 * R)
      );
      E = sqrt(
        (3 / 4) * a_3 ** 2 -
          R ** 2 -
          2 * a_2 -
          (4 * a_3 * a_2 - 8 * a_1 - a_3 ** 3) / (4 * R)
      );
    } else {
      D = sqrt((3 / 4) * a_3 ** 2 - 2 * a_2 + 2 * sqrt(y_1 ** 2 - 4 * a_0));
      E = sqrt((3 / 4) * a_3 ** 2 - 2 * a_2 - 2 * sqrt(y_1 ** 2 - 4 * a_0));
    }

    let r1 = -a_3 / 4 + R / 2 + D / 2;
    let r2 = -a_3 / 4 + R / 2 - D / 2;
    let r3 = -a_3 / 4 - R / 2 + E / 2;
    let r4 = -a_3 / 4 - R / 2 - E / 2;

    // Only return the real ones - complex ones will turn into NaN
    return [r1, r2, r3, r4].filter((r) => !Number.isNaN(r));
  }

  static transform(thing, mat = null) {
    // Vector
    if (thing instanceof p5.Vector) {
      let scaleFactor = 1;
      if (mat == null) {
        mat = drawingContext.getTransform();
        scaleFactor = pixelDensity();
      }

      let pt = mat.transformPoint(new DOMPoint(thing.x, thing.y));
      return createVector(pt.x / scaleFactor, pt.y / scaleFactor);
    }

    // Coord
    if (thing instanceof HitboxCoord) {
      return new HitboxCoord(Hitbox.transform(thing.p));
    }

    // Line
    if (thing instanceof HitboxLine) {
      return new HitboxLine(
        Hitbox.transform(thing.p1, mat),
        Hitbox.transform(thing.p2, mat)
      );
    }

    // Polygon
    if (thing instanceof HitboxPoly) {
      return new HitboxPoly(
        thing.ps.map((p) => Hitbox.transform(p, mat)),
        thing.includeFill
      );
    }

    // Ellipse
    if (thing instanceof HitboxEllipse) {
      return new HitboxEllipse(
        Hitbox.transform(thing.center, mat),
        Hitbox.transform(thing.dia1, mat),
        Hitbox.transform(thing.dia2, mat),
        thing.includeFill
      );
    }

    // Arc
    if (thing instanceof HitboxArc) {
      return new HitboxArc(
        Hitbox.transform(thing.ellipse, mat),
        Hitbox.transform(thing.startV, mat),
        Hitbox.transform(thing.stopV, mat),
        thing.mode
      );
    }
  }

  // Private static collision functions
  // ==================================

  // Ray casting algorithm - see https://stackoverflow.com/a/218081/2674563
  static coordInsidePoly(coord, poly) {
    return _hitbox_canvas
      .getContext("2d")
      .isPointInPath(poly.path, coord.p.x, coord.p.y);
  }

  static coordInsideEllipse(coord, ellipse) {
    // Geometric definition of an ellipse:
    // The set of all points for which the sum of a point's distance from the two foci is equal to some constant
    return (
      ellipse.f1.dist(coord.p) + ellipse.f2.dist(coord.p) <= ellipse.r + EPSILON
    );
  }

  static coordInsideArc(coord, arc) {
    // Firstly, the coordinate should definitely be within the entire ellipse
    if (!Hitbox.coordInsideEllipse(coord, arc.ellipse)) {
      return false;
    }

    // The coordinate's angle has to be between the start and stop angles
    let angle = Math.atan2(
      coord.p.y - arc.ellipse.center.y,
      coord.p.x - arc.ellipse.center.x
    );

    // Handling for chord-style arcs (CHORD or OPEN mode):
    // We simply check which side of the chord-line the point is on.
    if (arc.mode == p5.prototype.CHORD || arc.mode == p5.prototype.OPEN) {
      // Calculate arc length
      let start = arc.start - TWO_PI * Math.floor(arc.start / TWO_PI);
      let stop = arc.stop - TWO_PI * Math.floor(arc.stop / TWO_PI);
      let arcLength = Math.abs(start - stop);

      // If it's a small chord (arc length < PI), the point and ellipse center should be on opposite sides of the chord-line
      if (arcLength < PI) {
        return (
          arc.lines[0].sideOf(coord.p) !=
          arc.lines[0].sideOf(arc.ellipse.center)
        );
      }
      // If it's a big chord (arc length > PI), the point and ellipse center should be on the same side of the chord-line
      else if (arcLength > PI) {
        return (
          arc.lines[0].sideOf(coord.p) ==
          arc.lines[0].sideOf(arc.ellipse.center)
        );
      }
      // If it's a perfect half-chord (arc length = PI), it's identical to a non-chord arc, so fall through to the rest of handling
    }

    // Handling for pie-style arcs (null or PIE mode):
    // We simply check the angle
    return Hitbox.isAngleBetween(angle, arc.start, arc.stop);
  }

  static collideCoordCoord(coord1, coord2) {
    return coord1.p.dist(coord2.p) < EPSILON;
  }

  static collideCoordLine(coord, line) {
    return (
      Math.abs(
        coord.p.dist(line.p1) + coord.p.dist(line.p2) - line.p1.dist(line.p2)
      ) < EPSILON
    );
  }

  static collideCoordPoly(coord, poly) {
    return poly.lines.some((l) => Hitbox.collideCoordLine(coord, l));
  }

  static collideCoordEllipse(coord, ellipse) {
    // Geometric definition of an ellipse:
    // The set of all points for which the sum of a point's distance from the two foci is equal to some constant
    return (
      Math.abs(
        ellipse.f1.dist(coord.p) + ellipse.f2.dist(coord.p) - ellipse.r
      ) <= EPSILON
    );
  }

  static collideCoordArc(coord, arc) {
    // The coordinate's angle has to be between the start and stop angles
    let angle = Math.atan2(
      coord.p.y - arc.ellipse.center.y,
      coord.p.x - arc.ellipse.center.x
    );

    return (
      (Hitbox.isAngleBetween(angle, arc.start, arc.stop) &&
        Hitbox.collideCoordEllipse(coord, arc.ellipse)) ||
      arc.lines.some((l) => Hitbox.collideCoordLine(coord, l))
    );
  }

  static collideLineLine(line1, line2, returnCollisions) {
    let denom = line1.v.x * line2.v.y - line2.v.x * line1.v.y;

    if (denom == 0) {
      // Lines are parallel, so we just need to check if either of one line's endpoints are within the other
      if (
        Hitbox.collideCoordLine(new HitboxCoord(line1.p1), line2) ||
        Hitbox.collideCoordLine(new HitboxCoord(line1.p2), line2)
      ) {
        // Infinite collisions (unless we're just barely adjacent, but whatever)
        return returnCollisions ? [Infinity] : true;
      } else {
        return returnCollisions ? [] : false;
      }
    }

    let t1 =
      (line2.v.x * (line1.p1.y - line2.p1.y) +
        line2.v.y * (line2.p1.x - line1.p1.x)) /
      denom;
    let intersection = new HitboxCoord(
      p5.Vector.add(line1.p1, p5.Vector.mult(line1.v, t1))
    );

    let didCollide =
      Hitbox.collideCoordLine(intersection, line1) &&
      Hitbox.collideCoordLine(intersection, line2);

    if (returnCollisions) {
      return didCollide ? [intersection] : [];
    } else {
      return didCollide;
    }
  }

  static collideLinePoly(line, poly) {
    return poly.lines.some((l) => Hitbox.collideLineLine(l, line));
  }

  // If all you care about is whether an intersection exists,
  // set transformBack to false to save the computation time
  static getLineEllipseIntersection(line, ellipse, transformBack = true) {
    // Translate and rotate the line so we can put the ellipse flat at the origin, for cleaner math
    let mat = new DOMMatrix()
      .rotate((180 / PI) * -ellipse.angle)
      .translate(-ellipse.center.x, -ellipse.center.y);
    let tline = Hitbox.transform(line, mat);

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

      if (!transformBack) {
        return (t1 >= 0 && t1 <= tline.t) || (t2 >= 0 && t2 <= tline.t);
      }

      let intersections = [];
      for (const t of [t1, t2]) {
        if (t >= 0 && t <= tline.t) {
          intersections.push(
            Hitbox.transform(
              p5.Vector.add(tline.p1, p5.Vector.mult(tline.v, t)),
              mat.inverse()
            )
          );
        }
      }
      return intersections;
    }

    return transformBack ? [] : false;
  }

  static collideLineEllipse(line, ellipse) {
    return Hitbox.getLineEllipseIntersection(line, ellipse, false);
  }

  static collideLineArc(line, arc) {
    let intersections = Hitbox.getLineEllipseIntersection(line, arc.ellipse);
    return (
      intersections.some((p) =>
        Hitbox.collideCoordArc(new HitboxCoord(p), arc)
      ) || arc.lines.some((l) => Hitbox.collideLineLine(l, line))
    );
  }

  static collidePolyPoly(poly1, poly2) {
    for (const line1 of poly1.lines) {
      for (const line2 of poly2.lines) {
        if (Hitbox.collideLineLine(line1, line2)) {
          return true;
        }
      }
    }
    return false;
  }

  static collidePolyEllipse(poly, ellipse) {
    return poly.lines.some((l) => Hitbox.collideLineEllipse(l, ellipse));
  }

  static collidePolyArc(poly, arc) {
    return poly.lines.some((l) => Hitbox.collideLineArc(l, arc));
  }

  // If all you care about is whether an intersection exists,
  // set transformBack to false to save the computation time
  static getEllipseEllipseIntersection(
    ellipse1,
    ellipse2,
    transformBack = true
  ) {
    /*
    1. Move plane so ellipse1 is centered at the origin
    2. Rotate plane so ellipse1 is flat (major axis parallel to x axis)
    3. Transform plane so ellipse1 is a unit circle
    (order of operations is reversed due to how matrix math works)
    */
    let mat1 = new DOMMatrix()
      .scale(1 / ellipse1.a, 1 / ellipse1.b)
      .rotate((180 / PI) * -ellipse1.angle)
      .translate(-ellipse1.center.x, -ellipse1.center.y);
    ellipse1 = Hitbox.transform(ellipse1, mat1);
    ellipse2 = Hitbox.transform(ellipse2, mat1);

    /*
    4. Move plane so ellipse2 is centered at the origin
    5. Rotate plane so ellipse2 is flat (major axis parallel to x axis)
    */
    let mat2 = new DOMMatrix()
      .rotate((180 / PI) * -ellipse2.angle)
      .translate(-ellipse2.center.x, -ellipse2.center.y);
    ellipse1 = Hitbox.transform(ellipse1, mat2);
    ellipse2 = Hitbox.transform(ellipse2, mat2);

    // Rename important variables for easier use
    let m = ellipse1.center.x;
    let n = ellipse1.center.y;
    let w = ellipse2.a;
    let z = ellipse2.b;

    // Calculate coefficients of intersection quartic
    let z_0 = (1 - m ** 2 - n ** 2 - w ** 2) ** 2 - 4 * m ** 2 * w ** 2;
    let z_1 = 4 * n * (1 - m ** 2 - n ** 2 - w ** 2);
    let z_2 =
      2 * (m ** 2 + 3 * n ** 2 - 1) +
      ((2 * w ** 2) / z ** 2) * (m ** 2 - n ** 2 + z ** 2 - w ** 2 + 1);
    let z_3 = 4 * n * (w ** 2 / z ** 2 - 1);
    let z_4 = (w ** 2 / z ** 2 - 1) ** 2;

    // Solve the quartic for the y values of all intersections
    let ys = Hitbox.solveQuartic(z_4, z_3, z_2, z_1, z_0);

    // Filter identical solutions (multiple roots)
    let temp = [];
    outerloop: for (const y1 of ys) {
      for (const y2 of temp) {
        if (Math.abs(y1 - y2) < EPSILON) {
          continue outerloop;
        }
      }
      temp.push(y1);
    }
    ys = temp;

    // Recover the x values associated with each found y value
    let points = [];
    for (const y of ys) {
      // We have two possible options for X
      let x1 = m + Math.sqrt(1 - (y - n) ** 2);
      let x2 = m - Math.sqrt(1 - (y - n) ** 2);
      // The right one squared will be equal to this
      let xSq = w ** 2 - (w ** 2 / z ** 2) * y ** 2;

      // Pick the correct option for X
      let x;
      if (Math.abs(xSq) < EPSILON) {
        // Special case: an X value very near 0 can mess with the calculation due to impercision. Just set to 0
        x = 0;
      } else if (Math.abs(x1 ** 2 - xSq) < EPSILON) {
        x = x1;
      } else if (Math.abs(x2 ** 2 - xSq) < EPSILON) {
        x = x2;
      } else {
        // KNOWN ISSUE: for one-point intersections, we sometimes get incorrect solutions for y.
        // These are filtered here.
        continue;
      }

      if (!transformBack) {
        points.push(createVector(x, y));
      } else {
        // Transform the solution back into the original coordinate space by undoing the two transforms we did at the start
        let p = mat2
          .multiply(mat1)
          .inverse()
          .transformPoint(new DOMPoint(x, y));
        points.push(createVector(p.x, p.y));
      }
    }

    return points;
  }

  static collideEllipseEllipse(ellipse1, ellipse2) {
    // We can just check if any intersections exist without worrying about where they are
    return (
      Hitbox.getEllipseEllipseIntersection(ellipse1, ellipse2, false).length > 0
    );
  }

  static collideEllipseArc(ellipse, arc) {
    let intersections = Hitbox.getEllipseEllipseIntersection(
      ellipse,
      arc.ellipse
    );

    // At least one intersection has to be on the arc, or the arc's lines (if any) have to inersect with the ellipse
    return (
      intersections.some((p) =>
        Hitbox.collideCoordArc(new HitboxCoord(p), arc)
      ) || arc.lines.some((l) => Hitbox.collideLineEllipse(l, ellipse))
    );
  }

  static collideArcArc(arc1, arc2) {
    let intersections = Hitbox.getEllipseEllipseIntersection(
      arc1.ellipse,
      arc2.ellipse
    );

    // Either one intersection is on both arcs, a line hits an arc, or a line hits a line
    return (
      intersections.some(
        (p) =>
          Hitbox.collideCoordArc(new HitboxCoord(p), arc1) &&
          Hitbox.collideCoordArc(new HitboxCoord(p), arc2)
      ) ||
      arc1.lines.some((l) => Hitbox.collideLineArc(l, arc2)) ||
      arc2.lines.some((l) => Hitbox.collideLineArc(l, arc1)) ||
      arc1.lines.some((l1) =>
        arc2.lines.some((l2) => Hitbox.collideLineLine(l1, l2))
      )
    );
  }
}

// Capture mode overrides
Hitbox.prototype.currentCapturingHitbox = null;
let _wrapGraphicsFunctionForCapture = function (fName) {
  let f_base = p5.prototype[fName];
  p5.prototype[fName] = function (...args) {
    let curr = Hitbox.prototype.currentCapturingHitbox;
    if (curr != null) {
      curr[fName].apply(curr, [...args, DRAW_OFF]);
    }
    return f_base.apply(p5.instance, args);
  };
};
for (const fName of [
  "point",
  "line",
  "rect",
  "triangle",
  "quad",
  "ellipse",
  "arc",
  "square",
  "circle",
]) {
  _wrapGraphicsFunctionForCapture(fName);
}
