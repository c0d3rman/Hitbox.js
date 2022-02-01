# Hitbox.js

A library to provide simple 2D hitbox collision detection in p5.js programs.

Add to your p5js.org project by adding this line to index.html:

```html
<script src="https://cdn.jsdelivr.net/gh/c0d3rman/p5-libraries@master/Hitbox.min.js"></script>
```

Hitbox.js supports all the 2D primitives from p5.js: points, lines, ellipses, arcs, and polygons of all shapes.
To create a Hitbox object, simply use p5 transformations and functions as normal, and 'capture' them in a manner similar to `push()` and `pop()`

```javascript
h1 = new Hitbox();
h1.startCapture();
line(0, 0, 100, 100);
triangle(10, 10, 50, 70, 100, 40);
push();
rotate(0.5);
translate(100, 100);
ellipse(100, 100, 100, 200);
pop();
h1.stopCapture()
```

You can also create a hitbox directly without drawing anything:
```javascript
h2 = new Hitbox();
h2.line(0, 0, 100, 100);
h2.triangle(10, 10, 50, 70, 100, 40);
push();
rotate(0.5);
translate(100, 100);
h2.ellipse(100, 100, 100, 200);
pop();
```

Then check if a hitbox collides with another or with a point:

```javascript
if (h1.doesCollide(h2)) {
  // Do stuff
}
if (h1.doesCollidePoint(mouseX, mouseY)) {
  // Do stuff
}
```
