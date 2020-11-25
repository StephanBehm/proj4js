import adjust_lon from '../common/adjust_lon';
import HALF_PI from '../constants/values';
import proj4 from '../core';

/*
Original projection implementation:
https://github.com/OSGeo/PROJ/blob/46c47e9adf6376ae06afabe5d24a0016a05ced82/src/projections/ob_tran.cpp

Documentation:
https://proj.org/operations/projections/ob_tran.html

References/Formulas:
https://pubs.usgs.gov/pp/1395/report.pdf

Examples:
+proj=ob_tran +o_proj=moll +o_lat_p=45 +o_lon_p=-90
+proj=ob_tran +o_proj=moll +o_lat_p=45 +o_lon_p=-90 +lon_0=60
+proj=ob_tran +o_proj=moll +o_lat_p=45 +o_lon_p=-90 +lon_0=-90
*/
const projectionType = {
    not_set: 0,
    oblique: 1,
    transverse: 2 
};

export function init() {
    //1. Set parameters
    //Required
    this.o_proj = this.o_proj || "not_set" //Oblique projection - default to moll projection
    
    //New pole
    this.o_lat_p = this.o_lat_p || 0; //Latitude of the North pole of the unrotated source CRS, expressed in the rotated geographic CRS
    this.o_lon_p = this.o_lon_p || 0; //Longitude of the North pole of the unrotated source CRS, expressed in the rotated geographic CRS

    //Rotate about point
    this.o_alpha = this.o_alpha || 0; //Angle to rotate the projection with.
    this.o_lon_c = this.o_lon_c || 0; //Longitude of the point the projection will be rotated about.
    this.o_lat_c = this.o_lat_c || 0; //Latitude of the point the projection will be rotated about.

    //New “equator” points
    this.lon_1 = this.lon_1 || 0; //Longitude of first point.
    this.lat_1 = this.lat_1 || 0; //Latitude of first point.
    this.lon_2 = this.lon_2 || 0; //Longitude of second point.
    this.lat_2 = this.lat_2 || 0; //Latitude of second point.

    //Optional
    this.lon_0 = this.lon_0 || 0; //Longitude of projection center
    this.R = this.R; //Radius of the sphere given in meters. If used in conjunction with +ellps +R takes precedence.
    this.x_0 = this.x_0 || 0; //False easting
    this.y_0 = this.y_0 || 0; //False northing

    this.title = this.title || "General Oblique Transformation"; 

    this.TOL = 1e-10

    var Q; // Contains opaque/internal members (See L.11)
    Q.projectionType = projectionType.not_set;
    var phip; //double
    
    //Check inner projection
    if(this.o_proj = "not_set")
        throw new Error("No projection to rotate set: " + this.o_proj)

    //Generate PROJ string
    var innerProjection = "+proj=" + projectionName; 

    //Default input projection
    var defaultProjection = "+proj=longlat +datum=WGS84 +no_defs"

    //Create rotation
    if(this.o_alpha != 0) {
        var lamc = parseFloat(this.o_lon_c)
        var phic = parseFloat(this.o_lat_c)
        var alpha = parseFloat(this.o_alpha)

        if(Math.abs(Math.abs(phic) - HALF_PI) <= this.TOL) {
            throw new Error("Latitude is 0 or Alpha equals 90")
        }
        //L214
        Q.lamp = lamc + Math.atan2(-Math.cos(alpha), -Math.sin(alpha) * Math.sin(phic));
        phip = Math.asin(Math.cos(phic) * Math.sin(alpha));
    }
    else if(o_lat_p != 0)  { //New pole specified
        //L216
        Q.lamp = parseFloat(this.o_lon_p);
        phip = parseFloat(this.o_lat_p);
    }
    else { //New equator points specified
        //L219
        var lam1, lam2, phi1, phi2, con;

        lam1 = parseFloat(this.lon_1);
        phi1 = parseFloat(this.lat_1);
        lam2 = parseFloat(this.lon_2);
        phi2 = parseFloat(this.lat_2);

        con = Math.abs(phi1)

        if(Math.abs(phi1 - phi2) <= this.TOL || con <= TOL || 
            Math.abs(con - HALF_PI) <= TOL || Math.abs(Math.abs(phi2) - HALF_PI) <= TOL)
                throw new Error("Lat1 or Lat2 are 0 or 90")

        Q.lamp = Math.atan2(Math.cos(phi1) * Math.sin(phi2) * Math.cos(lam1) 
                - Math.sin(phi1) * Math.cos(phi2) * Math.cos(lam2), 
                Math.sin(phi1) * Math.cos(phi2) * Math.sin(lam2) 
                - Math.cos(phi1) * Math.sin(phi2) * Math.sin(lam1));
        phip = Math.atan(-Math.cos(Q.lamp - lam1) / Math.tan(phi1));        
    }

    if(Math.abs(phip) > this.TOL) {
        console.log("Oblique projection")

        Q.cphip = Math.cos(phip);
        Q.sphip = Math.sin(phip)

        Q.projectionType = projectionType.oblique;        
    } else {
        console.log("Transverse projection")
        Q.projectionType = projectionType.transverse;
    }
}

// forward equations--mapping (lat,long) to (x,y)
// oblique - true poles of earth lie on the equator of the basic projection,
//  and the poles of the projection lie on the equator of the earth
// transverse -
//   
// -----------------------------------------------------------------
export function forward(p) {
    //var Q = pj_opaque;
    
    var lam = p.x; //Lambda
    var phi = p.y; //Phi

    //Oblique or transverse
    switch(Q.projectionType) {
        case projectionType.not_set: 
            throw new Error("ProjectionType not set on Forward-Transformation")

        case projectionType.oblique: {
            var coslam = Math.cos(lam)
            var sinphi = Math.sin(phi)
            var cosphi = Math.cos(phi)

            /* Formula (5-8b) of Snyder's "Map projections: a working manual" */
            lam = adjust_lon(Math.atan2(cosphi * Math.sin(lam), Q.sphip * cosphi * coslam + Q.cphip * sinphi) + Q.lamp);
            /* Formula (5-7) */
            phi = Math.asin(Q.sphip * sinphi - Q.cphip * cosphi * coslam);  

            p.x = lam
            p.y = phi

            console.log("Before inner projection: " + p)
            p = proj4(defaultProjection,innerProjection).forward(p.toPoint())
            console.log("After inner projection: " + p)


            return p;          
        }

        case projectionType.transverse: {
            var cosphi = Math.cos(phi)
            var coslam = Math.cos(lam)

            lam = adjust_lon(Math.atan2(cosphi * Math.sin(lam), Math.sin(phi)) + Q.lamp);
            phi = Math.asin(-cosphi * coslam);

            p.x = lam
            p.y = phi

            console.log("Before inner projection: " + p)
            p = proj4(defaultProjection,innerProjection).forward(p.toPoint())
            console.log("After inner projection: " + p)

            return p;  
        }

        default:
            throw new Error("Unknown ProjectionType")
    }     
}

// inverse equations--mapping (x,y) to (lat,long) 

export function inverse(p) {
    
    switch(Q.projectionType) {
        case projectionType.not_set: 
            throw new Error("ProjectionType not set on Forward-Transformation")

        case projectionType.oblique: {
            //TODO::Implement
            var coslam = Math.cos(lam)
            var sinphi = Math.sin(phi)
            var cosphi = Math.cos(phi)

            //Run inverse projection
            console.log("Before inner projection: " + p)
            p = proj4(innerProjection,defaultProjection).inverse(p.toPoint())
            console.log("After inner projection: " + p)

            var lam = p.x; //Lambda
            var phi = p.y; //Phi

            //L63
            if(lam < Number.MAX_VALUE) {
                lam -= Q.lamp;
                coslam = Math.cos(lam);
                sinphi = Math.sin(phi);
                cosphi = Math.cos(phi)
                /* Formula (5-9) */
                phi = Math.asin(Q.sphip * sinphi + Q.cphip * cosphi * coslam);
                /* Formula (5-10b) */
                lam = Math.atan2(cosphi * Math.sin(lam),Q.sphip * cosphi * coslam - Q.cphip * sinphi)
            }
            p.x = lam
            p.y = phi

            return p;
        }

        case projectionType.transverse: {
            //TODO:

        }

        default:
            throw new Error("Unknown ProjectionType")
    }   
}

export var names = ["General Oblique Transformation","General_Oblique_Transformation","obtran"];
export default {
  init: init,
  forward: forward,
  inverse: inverse,
  names: names
};