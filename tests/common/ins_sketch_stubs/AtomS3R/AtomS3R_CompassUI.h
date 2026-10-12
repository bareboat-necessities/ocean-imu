#pragma once
// Host stand-in for the sprite compass renderer.
class CompassUI {
 public:
  void begin() {}
  bool ok() const { return false; }
  template <class... T> void draw(T...) {}
};
