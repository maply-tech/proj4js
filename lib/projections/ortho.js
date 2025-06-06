import adjust_lon from '../common/adjust_lon';
import { EPSLN } from '../constants/values';

const M_PI = Math.PI;
const M_TWOPI = 2. * M_PI;
const M_PI_2 = M_PI / 2.;

const adjlon = (l) => {
  if (Math.abs(l) < M_PI + 1e-12) return l;

  l += M_PI;
  l -= M_TWOPI * Math.floor(l / M_TWOPI);
  l -= M_PI;

  return l;
}

export function init() {
  // double temp;      /* temporary variable    */

  /* Place parameters in static storage for common use
    ------------------------------------------------- */
  this.sinph0 = Math.sin(this.lat0);
  this.cosph0 = Math.cos(this.lat0);
  this.nu0 = 1. / Math.sqrt(1. - this.es * this.lat0 * this.lat0);
}

/* Orthographic forward equations--mapping lat,long to x,y
    --------------------------------------------------- */
export function forward(p) {
  var sinphi, cosphi; /* sin and cos value        */
  var dlon;
  var sinlam, coslam; /* cos of longitude        */
  var x, y;
  var lon = p.x;
  var lat = p.y;
  /* Forward equations
      -----------------*/
  dlon = adjust_lon(lon - this.long0);

  sinphi = Math.sin(lat);
  cosphi = Math.cos(lat);
  sinlam = Math.sin(dlon);
  coslam = Math.cos(dlon);

  var nu = 1. / Math.sqrt(1. - this.es * sinphi * sinphi);
  x = nu * cosphi * sinlam;
  y = nu * (sinphi * this.cosph0 - cosphi * this.sinph0 * coslam) +
    this.es * (this.nu0 * this.sinph0 - nu * sinphi) * this.cosph0;

  // g = this.sinph0 * sinphi + this.cosph0 * cosphi * coslam;
  // if ((g > 0) || (Math.abs(g) <= EPSLN)) {
  //   x = cosphi * Math.sin(dlon);
  //   y = (this.cosph0 * sinphi - this.sinph0 * cosphi * coslam);
  // }
  p.x = x;
  p.y = y;

  var alpha = this.alpha || 0;
  var cosalpha = Math.cos(alpha);
  var sinalpha = Math.sin(alpha);
  p.x = (x * cosalpha - y * sinalpha) * this.a + this.x0;
  p.y = (x * sinalpha + y * cosalpha) * this.a + this.y0;

  return p;
}

export function inverse(p) {
  var rh;
  var phi, lam;

  /* Inverse equations
      -----------------*/

  p.x -= this.x0;
  p.y -= this.y0;

  var x = p.x / this.a / this.k0;
  var y = p.y / this.a / this.k0;

  var alpha = this.alpha || 0;
  var cosalpha = Math.cos(alpha);
  var sinalpha = Math.sin(alpha);

  let xy = {
    x: (x * cosalpha + y * sinalpha),
    y: (-x * sinalpha + y * cosalpha),
  };

  const c = { x: xy.x, y: xy.y };

  rh = Math.sqrt(xy.x * xy.x + xy.y * xy.y);

  if (Math.abs(rh) <= EPSLN) {
    p.x = this.long0;
    p.y = this.lat0;
    return p;
  }

  const sinc = rh > 1. ? 1. : rh;
  const cosc = Math.sqrt(1. - sinc * sinc);

  phi = cosc * this.sinph0 + xy.y * sinc * this.cosph0 / rh;
  xy.y = (cosc - this.sinph0 * phi) * rh;
  xy.x *= sinc * this.cosph0;
  phi = Math.asin(phi);
  lam = Math.atan2(xy.x, xy.y);


  for (let i = 0; i < 20; i++) {
    var cosphi = Math.cos(phi);
    var sinphi = Math.sin(phi);
    var coslam = Math.cos(lam);
    var sinlam = Math.sin(lam);
    var one_minus_es_sinphi2 = 1. - this.es * sinphi * sinphi;
    var nu = 1. / Math.sqrt(one_minus_es_sinphi2);

    var x = nu * cosphi * sinlam;
    var y = nu * (sinphi * this.cosph0 - cosphi * this.sinph0 * coslam) +
      this.es * (this.nu0 * this.sinph0 - nu * sinphi) * this.cosph0;

    var rho = (1. - this.es) * nu / one_minus_es_sinphi2;
    var J11 = -rho * sinphi * sinlam;
    var J12 = nu * cosphi * coslam;
    var J21 = rho * (cosphi * this.cosph0 + sinphi * this.sinph0 * coslam);
    var J22 = nu * this.cosph0 * cosphi * sinlam;
    var D = J11 * J22 - J12 * J21;
    var dx = c.x - x;
    var dy = c.y - y;
    var dphi = (J22 * dx - J12 * dy) / D;
    var dlam = (-J21 * dx + J11 * dy) / D;
    if (phi > M_PI_2) {
      phi = M_PI_2 - (phi - M_PI_2);
      lam = adjlon(lam + M_PI);
    } else if (phi < -M_PI_2) {
      phi = -M_PI_2 + (-M_PI_2 - phi);
      lam = adjlon(lam + M_PI);
    }
    phi += dphi;
    lam += dlam;
    if (Math.abs(dphi) < 1e-12 && Math.abs(dlam) < 1e-12) break;
  }

  p.x = lam + this.long0;
  p.y = phi;

  return p;
}

export var names = ['ortho'];
export default {
  init: init,
  forward: forward,
  inverse: inverse,
  names: names
};
