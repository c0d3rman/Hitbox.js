# p5-libraries

A collection of libraries to be used with p5js.
It is recommended to use jsdelivr to load these.

## Hitbox.js

Add to your p5js.org project by adding this line to index.html:

```html
<script src="https://cdn.jsdelivr.net/gh/c0d3rman/p5-libraries@master/Hitbox.min.js"></script>
```

A library to provide simple 2D hitbox collision detection. Supports points, lines, polygons of all shapes, and ellipses.
Simply create a Hitbox object, using p5 transformations and functions as normal:

```javascript
h1 = new Hitbox();
h1.line(0, 0, 100, 100);
h1.triangle(10, 10, 50, 70, 100, 40);
push();
rotate(0.5);
translate(100, 100);
h1.ellipse(100, 100, 100, 200);
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
