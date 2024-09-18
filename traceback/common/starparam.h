// $Id: starparam.h 143 2022-07-28 10:52:27Z ifg $

#ifndef STARPARAM_H
#define STARPARAM_H

#include <QStringList>
#include <QFile>
#include <QTextStream>
#include <QSet>
#include "Matrix.h"
#include "rng.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#define CRITHOOGERWERF "Hoogerwerf"
#define CRITTETZLAFF "Tetzlaff"
#define CRITRALPH1 "Ralph1"
#define CRITVARLERI1 "Valeri1"

/*!
 * \mainpage
 * Title page.
 *
 * It should contain a few words about the
 * coordinate systems used, and the transformations
 * between them.
 *
 * Some day it will.
 */

typedef double mytype;

/*!
 * \brief The GalaxyParams struct
 * contains parameters to describe the galactic potential
 * of the milky way.
 *
 * The values are read from the configuration file.
 * \sa StarGalacticParam::eqnOfMotion
 */
struct GalaxyParams {
	double G,
	       Mdisk,
	       Msphere,
	       vhalo,
	       a, b, c, d;
};

/*!
 * \brief The SimulParams struct
 * contains the parameters for the simulation
 *
 * The values are read from the configuration file
 */
struct SimulParams {
	int     Threads,	///< Number of threads to run in parallel (0 = auto)
	        Orbits,		///< Number of orbits to calculate
	        Steps,		///< Number of steps for each orbit
	        Width;		///< Step width for each orbit (in years)
	mytype  StepSize,	///< Initial step size for integration
	        Limit;		///< Upper limit of the minimum separation that is considered a success (in parsec)
	bool    RemoveFile;	///< Remove %output file if no successful orbits were found
	QString Star1,		///< FileName\#LineNr for first star
	        Star2,		///< FileName\#LineNr for second star
	        Assoc;		///< FileName\#LineNr for star association (0 = no association)
	QString Type,		///< Orbit type to use ('Potential', 'Epicycle', or 'Linear')
	        Criterion,	///< Which orbits are considered successful (see comments in config file)
	        InFile,		///< Input file with stellar co-ordinates
	        OutFile;	///< %Output file for results, if omitted or set to "auto", basename of configuration file with extension '.out' will be used
};
/* FIXME: How can QStrings in a structure work? */

/**
 * @brief Class to read and hold the star parameters from a file
 *
 * The StarParam class reads and holds values from an input file.
 * This file must have one line per star. The first entry of
 * each line designates the type on which the number and
 * meaning of the following values depend.
 * <b>The order is important.</b>
 *
 * \sa #StarParam::CoordType
 * \sa #StarParam::DistType
 *
 * Supported types are (so far):
 *
 * \li \b 1 - Equatorial coordinates; at least 17 values having the following meaning:
 *
 * -# Type of values (= \b 1)
 * -# Right ascension [hh]
 * -# Right ascension [mm]
 * -# Right ascension [ss.sss]
 * -# Declination [°]
 * -# Declination [']
 * -# Declination ["]
 * -# Parallax [mas]
 * -# Parallax error [mas]
 * -# Radial velocity [km/s]
 * -# Radial velocity error [km/s]
 * -# The distribution used to vary radial velocity (#DistType)
 * -# Proper motion in RA, corrected for declination
 *    @f$(\mu_{\alpha}^\star = \mu_{\alpha} \cos \delta)@f$ [mas/yr]
 * -# Proper motion error in RA [mas/yr]
 * -# Proper motion in Dec [mas/yr]
 * -# Proper motion error in Dec [mas/yr]
 * -# Identifier string or HIP# (no spaces!)
 * -# Optional: remarks (may contain spaces)
 *
 * \li \b 2 - GalacticUVW coordinates; at least 12 values having the following meaning:
 *
 * -# Type of values (= \b 2)
 * -# Galactic longitude [°]
 * -# Galactic latitude [°]
 * -# Parallax [mas]
 * -# Parallax error [mas]
 * -# Velocity U [km/s]
 * -# Error of U [km/s]
 * -# Velocity V [km/s]
 * -# Error of V [km/s]
 * -# Velocity W [kms]
 * -# Error of W [km/s]
 * -# Identifier string or HIP# (no spaces!)
 * -# Optional: remarks (may contain spaces)
 *
 * \li \b 3 - GalacticPM coordinates; at least 13 values having the following meaning:
 *
 * -# Type of values (= \b 3)
 * -# Galactic longitude [°]
 * -# Galactic latitude [°]
 * -# Parallax [mas]
 * -# Parallax error [mas]
 * -# Radial velocity [km/s]
 * -# Radial velocity error [km/s]
 * -# The distribution used to vary radial velocity (#DistType)
 * -# Proper motion in galactic longitude, corrected for latitude
 *    @f$(\mu_{\ell}^\star = \mu_{\ell} \cos b)@f$ [mas/yr]
 * -# Proper motion error in galactic longitude [mas/yr]
 * -# Proper motion in galactic latitude [mas/yr]
 * -# Proper motion error in galactic latitude [mas/yr]
 * -# Identifier string or HIP# (no spaces!)
 * -# Optional: remarks (may contain spaces)
 *
 * \li \b 4 - Equatorial coordinates with proper motions in ecliptic ccordinates;
 * at least 17 values having the following meaning:
 *
 * -# Type of values (= \b 4)
 * -# Right ascension [hh]
 * -# Right ascension [mm]
 * -# Right ascension [ss.sss]
 * -# Declination [°]
 * -# Declination [']
 * -# Declination ["]
 * -# Parallax [mas]
 * -# Parallax error [mas]
 * -# Radial velocity [km/s]
 * -# Radial velocity error [km/s]
 * -# The distribution used to vary radial velocity (#DistType)
 * -# Proper motion in ecliptic longitute, corrected for latitude
 *    @f$(\mu_{\lambda}^\star = \mu_{\lambda} \cos \beta)@f$ [mas/yr]
 * -# Proper motion error in ecliptic longitude [mas/yr]
 * -# Proper motion in ecliptic latitude [mas/yr]
 * -# Proper motion error in ecliptic latitude [mas/yr]
 * -# Identifier string or HIP# (no spaces!)
 * -# Optional: remarks (may contain spaces)
 *
 * \li \b 101 - Equatorial coordinates for association;
 * at least 18 values having the following meaning:
 *
 * -# Type of values (= \b 101)
 * -# Right ascension [hh]
 * -# Right ascension [mm]
 * -# Right ascension [ss.sss]
 * -# Declination [°]
 * -# Declination [']
 * -# Declination ["]
 * -# Parallax [mas]
 * -# Parallax error [mas]
 * -# Radial velocity [km/s]
 * -# Radial velocity error [km/s]
 * -# The distribution used to vary radial velocity (#DistType)
 * -# Proper motion in RA, corrected for declination
 *    @f$(\mu_{\alpha}^\star = \mu_{\alpha} \cos \delta)@f$ [mas/yr]
 * -# Proper motion error in RA [mas/yr]
 * -# Proper motion in Dec [mas/yr]
 * -# Proper motion error in Dec [mas/yr]
 * -# Radius of association [pc]
 * -# Identifier string or HIP# (no spaces!)
 * -# Optional: remarks (may contain spaces)
 *
 * \li \b 102 - GalacticUVW coordinates for associations;
 * at least 13 values having the following meaning:
 *
 * -# Type of values (= \b 102)
 * -# Galactic longitude [°]
 * -# Galactic latitude [°]
 * -# Parallax [mas]
 * -# Parallax error [mas]
 * -# Velocity U [km/s]
 * -# Error of U [km/s]
 * -# Velocity V [km/s]
 * -# Error of V [km/s]
 * -# Velocity W [kms]
 * -# Error of W [km/s]
 * -# Radius of association [pc]
 * -# Identifier string or HIP# (no spaces!)
 * -# Optional: remarks (may contain spaces)
 *
 * \li \b 103 - GalacticPM coordinates for associations;
 * at least 14 values having the following meaning:
 *
 * -# Type of values (= \b 103)
 * -# Galactic longitude [°]
 * -# Galactic latitude [°]
 * -# Parallax [mas]
 * -# Parallax error [mas]
 * -# Radial velocity [km/s]
 * -# Radial velocity error [km/s]
 * -# The distribution used to vary radial velocity {#DistType}
 * -# Proper motion in galactic longitude, corrected for latitude
 *    @f$(\mu_{\ell}^\star = \mu_{\ell} \cos b)@f$ [mas/yr]
 * -# Proper motion error in galactic longitude [mas/yr]
 * -# Proper motion in galactic latitude [mas/yr]
 * -# Proper motion error in galactic latitude [mas/yr]
 * -# Radius of association [pc]
 * -# Identifier string or HIP# (no spaces!)
 * -# Optional: remarks (may contain spaces)
 *
 * \li \b 104 - Equatorial coordinates with proper motions in ecliptic coordiantes for association;
 * at least 18 values having the following meaning:
 *
 * -# Type of values (= \b 104)
 * -# Right ascension [hh]
 * -# Right ascension [mm]
 * -# Right ascension [ss.sss]
 * -# Declination [°]
 * -# Declination [']
 * -# Declination ["]
 * -# Parallax [mas]
 * -# Parallax error [mas]
 * -# Radial velocity [km/s]
 * -# Radial velocity error [km/s]
 * -# The distribution used to vary radial velocity (#DistType)
 * -# Proper motion in ecliptic longitude, corrected for latitude
 *    @f$(\mu_{\lambda}^\star = \mu_{\lambda} \cos \beta)@f$ [mas/yr]
 * -# Proper motion error in ecpliptic longitude [mas/yr]
 * -# Proper motion in ecliptic latitude [mas/yr]
 * -# Proper motion error in ecliptic latitude [mas/yr]
 * -# Radius of association [pc]
 * -# Identifier string or HIP# (no spaces!)
 * -# Optional: remarks (may contain spaces)
 *
 * The values are converted as necessary.
 */
class StarParam
{
public:
	/*!
	 * \brief The CoordType enum
	 *
	 * This enumerator describes the supported input formats i.e.,
	 * the meaning of the values in a single line of the input file
	 * #fromFile.
	 */
	enum CoordType {
		Equatorial = 1,				///< 1, @f$\alpha(3), \delta(3), \pi, \sigma_\pi, v_{rad}, \sigma_{v_{rad}}@f$,
									///  DistType, @f$\mu_{\alpha}^*, \sigma_{\mu_{\alpha}^*},
									///  \mu_{\delta}, \sigma_{\mu_{\delta}}@f$, Id, Remarks
		GalacticUVW = 2,			///< 2, @f$\ell, b, \pi, \sigma_{\pi}, U, \sigma_U,
									///  V, \sigma_V, W, \sigma_W@f$ Id, Remarks
		GalacticPM = 3,				///< 3, @f$\ell, b, \pi, \sigma_{\pi}, v_{rad}, \sigma_{v_{rad}}@f$,
									///  DistType, @f$\mu_{\ell}^*, \sigma_{\mu_{\ell}^*},
									///  \mu_{b}, \sigma_{\mu_{b}}@f$, Id, Remarks
		EquatorialEPM = 4,			///< 1, @f$\alpha(3), \delta(3), \pi, \sigma_\pi, v_{rad}, \sigma_{v_{rad}}@f$,
									///  DistType, @f$\mu_{\lambda}^*, \sigma_{\mu_{\lambda}^*},
									///  \mu_{\beta}, \sigma_{\mu_{\beta}}@f$, Id, Remarks
		EquatorialAssoc = 101,		///< 101, @f$\alpha(3), \delta(3), \pi, \sigma_\pi, v_{rad}, \sigma_{v_{rad}}@f$,
									///  DistType, @f$\mu_{\alpha}^*, \sigma_{\mu_{\alpha}^*},
									///  \mu_{\delta}, \sigma_{\mu_{\delta}}, R@f$, Id, Remarks
		GalacticUVWAssoc = 102,		///< 102, @f$\ell, b, \pi, \sigma_{\pi}, U, \sigma_U,
									///  V, \sigma_V, W, \sigma_W, R@f$, Id, Remarks
		GalacticPMAssoc = 103,		///< 103, @f$\ell, b, \pi, \sigma_{\pi}, v_{rad}, \sigma_{v_{rad}}@f$,
									///  DistType, @f$\mu_{\ell}^*, \sigma_{\mu_{\ell}^*},
									///  \mu_{b}, \sigma_{\mu_{b}}, R@f$, Id, Remarks
		EquatorialEPMAssoc = 104	///< 101, @f$\alpha(3), \delta(3), \pi, \sigma_\pi, v_{rad}, \sigma_{v_{rad}}@f$,
									///  DistType, @f$\mu_{\lambda}^*, \sigma_{\mu_{\lambda}^*},
									///  \mu_{\beta}, \sigma_{\mu_{\beta}}, R@f$, Id, Remarks
	};

	/*!
	 * \brief The DistType enum
	 *
	 * This enumeration describes the statistical method
	 * that is used to vary the radial velocity of a star/association
	 * during the simulation.
	 *
	 * \sa
	 * <a href="http://adsabs.harvard.edu/abs/2005MNRAS.360..974H" target="_new">
	 * G. Hobbs et al. (2005): A statistical study of 233 pulsar proper motions.
	 * \a MNRAS \b 360, 974-992
	 * </a>
	 */
	enum DistType {
		Gaussian = 0,		///< Use observed mean value and standard deviation for radial velocity (default)
		Maxwellian = 1,		///< Assume a Maxwellian distribution for the space velocity (see Hobbs et al. 2005)
							///  and subtract transversal velocity
		Uniform = 2			///< Assume a Uniform distribution -1000 ... +1000 km/s for the radial velocity
	};

private:
	CoordType coordType;		///< Type of values read in
	union {
	  mytype  RA,				///< Right ascension [°]
	          l;				///< Galactic longitude[°]
	};
	union {
	  mytype  Dec,				///< Declination [°]
	          b;				///< Galactic latitude [°]
	};
	// The values in the input file must have the units listed here
	mytype    para,				///< Parallax [mas], used for calculation
	          para0,			///< Measured parallax [mas], as read from file
	          p_err;			///< Parallax error [mas]
	union {
	  mytype  RV,				///< Radial velocity [km/s], used for calculation
			  U;				///< Velocity in X direction [km/s]
	};
	union {
	  mytype  RV0,				///< Measured radial velocity [km/s], as read from file
	          U0;				///< Original velocity in X direction [km/s], as read from file
	};
	union {
	  mytype  RV_err,			///< Radial velocity error [km/s]
	          U_err;			///< Error of U [km/s]
	};
	union {
	  mytype  PMa,				///< Proper motion in RA, corrected for declination [mas/yr], used for calculation
	          PMl,				///< Proper motion in galactic longitude, corrected for latitude @f$(\mu_\ell \cos b)@f$ [mas/yr]
	          PMlambda,			///< Proper motion in ecliptic longitude, corrected for latitude @f$(\mu_\lambda \cos \beta)@f$ [mas/yr]
	          V;				///< Velocity in Y direction [km/s]
	};
	union {
	  mytype  PMa0,				///< Measured proper motion in RA, corrected for declination [mas/yr], as read from file
	          PMl0,				///< Measured proper motion in gal. longitude, corrected for latitude [mas/yr], as read from file
	          PMlambda0,		///< Measured proper motion in ecl. longitude, corrected for latitude [mas/yr], as read from file
	          V0;				///< Original velocity in Y direction [km/s], as read from file
	};
	union {
	  mytype  PMa_err,			///< Make an educated guess [mas/yr]
	          PMl_err,			///< Error of PMl [mas/yr]
	          PMlambda_err,		///< Error of PMlambda [mas/yr]
	          V_err;			///< Error of V [km/s]
	};
	union {
	  mytype  PMd,				///< Proper motion in Dec [mas/yr], used for calulation
	          PMb,				///< Proper motion in galactic latitude [mas/yr], used for calulation
	          PMbeta,			///< Proper motion in ecliptic latitude [mas/yr], used for calulation
	          W;				///< Velocity in Z direction [km/s]
	};
	union {
	  mytype  PMd0,				///< Measured proper motion in Dec [mas/yr], as read from file
	          PMb0,				///< Measured proper motion in gal. latitude [mas/yr], as read from file
	          PMbeta0,			///< Measured proper motion in ecl. latitude [mas/yr], as read from file
	          W0;				///< Original velocity in Z direction [km/s], as read from file
	};
	union {
	  mytype  PMd_err,			///< [mas/yr]
	          PMb_err,			///< [mas/yr]
	          PMbeta_err,		///< [mas/yr]
	          W_err;			///< Error of W [km/s]
	};
	mytype    radius;			///< Radius of association [pc]
	QString   HIP,				///< HIP number
	          rem;				///< Comment, remarks
	QString   fromFile;			///< File where the values were read from
	int       lineNumber;		///< Line# in @a #fromFile
	DistType  distType;			///< What distribution to use to simulate the radial velocity (for neutron stars with unknown RV)

public:
	/*!
	 * \brief Values prone to statistical errors
	 */
	enum Uncertainty {
		Parallax,			///< Parallax (distance) of the star
		RadialVelocity,		///< Velocity of the star in viewing direction
		PMAlpha,			///< Proper motion in right ascension, corrected for declination
		PMDelta,			///< Proper motion in declination
		PMLong,				///< Proper motion in galactic longitude, corrected for latitude
		PMLat,				///< Proper motion in galactic latitude
		PMELong,			///< Proper motion in ecliptic longitude, corrected for latitude
		PMELat,				///< Proper motion in ecliptic latitude
		VelU,				///< Velocity in X direction
		VelV,				///< Velocity in Y direction
		VelW				///< Velocity in Z direction
	};
	explicit StarParam();
	int readValuesFromFile(QFile &file, int lineNumber);
	mytype varyValue(Uncertainty what, Rng &rng);
	bool varyAllValues(Rng &rng);

	/// @name Return individual components
	/// @{
	mytype    getRA() const            {return RA;}
	mytype    getDec() const           {return Dec;}
	mytype    getParallax() const      {return para;}
	mytype    getParallax0() const     {return para0;}
	mytype    getParaerr() const       {return p_err;}
	mytype    getRV() const            {return RV;}
	mytype    getRV0() const           {return RV0;}
	mytype    getRVerr() const         {return RV_err;}
	mytype    getPMa() const           {return PMa;}
	mytype    getPMa0() const          {return PMa0;}
	mytype    getPMaerr() const        {return PMa_err;}
	mytype    getPMd() const           {return PMd;}
	mytype    getPMd0() const          {return PMd0;}
	mytype    getPMderr() const        {return PMd_err;}
	mytype    getPMl() const           {return PMl;}
	mytype    getPMl0() const          {return PMl0;}
	mytype    getPMlerr() const        {return PMl_err;}
	mytype    getPMb() const           {return PMb;}
	mytype    getPMb0() const          {return PMb0;}
	mytype    getPMberr() const        {return PMb_err;}
	mytype    getPMlambda() const      {return PMlambda;}
	mytype    getPMlambda0() const     {return PMlambda0;}
	mytype    getPMlambdaerr() const   {return PMlambda_err;}
	mytype    getPMbeta() const        {return PMbeta;}
	mytype    getPMbeta0() const       {return PMbeta0;}
	mytype    getPMbetaerr() const     {return PMbeta_err;}
	QString   getHIP() const           {return HIP;}
	QString   getRem() const           {return rem;}
	QString   getInputFileName() const {return fromFile;}
	int       getLineNumber() const    {return lineNumber;}
	CoordType getCoordType() const     {return coordType;}
	DistType  getDistType() const      {return distType;}
	mytype    getLongitude() const     {return l;}
	mytype    getLatitude() const      {return b;}
	mytype    getU() const             {return U;}
	mytype    getV() const             {return V;}
	mytype    getW() const             {return W;}
	mytype    getUerr() const          {return U_err;}
	mytype    getVerr() const          {return V_err;}
	mytype    getWerr() const          {return W_err;}
	mytype    getRadius() const        {return radius;}
	bool      isAssoc() const          {return (coordType>100);}
	/// @}

	void printValues(QTextStream &out, bool prependHash) const;
	static mytype  deg2rad(mytype degrees);
	static mytype  rad2deg(mytype radian);
	static QString deg2dms(mytype degrees);
	static QString deg2hms(mytype degrees);
	static mytype  dms2deg(Vector<mytype> dms, bool isNegative);
	static mytype  hms2deg(Vector<mytype> hms, bool isNegative);
};

/**
 * @brief The StarGalacticParam class
 *
 * The StarGalacticParam class will hold all values to describe
 * the location of the star in the galactic coordinate system @a l, @a b,
 * the cartesian system @a X, @a Y, @a Z, and the space velocities @a U, @a V, @a W.
 */
class StarGalacticParam
{
public:
	/**
	 * @brief The method used for calculating the orbit
	 */
	enum OrbitType {
		Epicycle,		///< Use an analytical approximation (see StarGalacticParam::epicyclicOrbit()).
						///< Coordinates are relative to the Sun with the \a X axis \b always pointing to the
						///< galactic center, i.e. the coordinate system is rotating around the galaxy
		Numeric,		///< Numerically integrate the equation of motion (see StarGalacticParam::numericOrbit()).
						///< Coordinates are relative to the Sun at time 0
		Linear			///< Assume linear motion through space determined by the space velocity
						///< (see StarGalacticParam::linearOrbit()).
						///< This is not really an orbit at all but might be sufficient
						///< for short time intervals. Uses the same coordinate system
						///< as #Numeric orbit.
	};

private:
	/// @name References to Cartesian coordinates [pc]
	/// @{
	mytype& x() { return cartesian[0]; }		///< Positive towards the Galactic center
	mytype& y() { return cartesian[1]; }		///< Positive in Galactic rotation
	mytype& z() { return cartesian[2]; }		///< Positive towards North Galactic Pole

	/// @}

	/// @name References to space-velocity components [pc/yr]
	/// @{
	mytype& u() { return cartesian[3]; }		///< In \a X direction
	mytype& v() { return cartesian[4]; }		///< In \a Y direction
	mytype& w() { return cartesian[5]; }		///< In \a Z direction

	/// @}

	/// @name Position and velocity in cartesian coordinates in [pc] and [pc/yr], resp.
	/// @{
	mytype cartesian[6];	///< References to the cartesian coordinates and space velocities
							///< are available by x(), y(), z(), u(), v(), w().
							///< The coordinates can be relative to an inertial system
							///< or a rotating system, depending on
							///< the orbit calculation method used.
							///< \sa OrbitType

	/// @}

	/// @name Time for coordinates
	/// @{
	mytype time;		///< [yr]

	/// @}

	/// @name Starting values for orbits (might have undergone Monte Carlo variation)
	/// @{
	mytype x0, y0, z0,
		   u0, v0, w0;
	/// @}

	/// @name Uncertainties of cartesian coordinates
	/// @{
	mytype sigmaX0, sigmaY0, sigmaZ0,
		   sigmaX, sigmaY, sigmaZ;
	/// @}

	/// @name Uncertainties of velocities
	/// @{
	mytype sigmaU0, sigmaV0, sigmaW0,
		   sigmaU, sigmaV, sigmaW;
	/// @}

	/// @name Transformation Matrix
	/// @{
	const Matrix<mytype> T;
	/// @}

	/// @name For solving the equation of motion numerically
	/// @{
	gsl_odeiv2_system odeSys;		///< The differential equation system (see eqnOfMotion() and jacEqnOfMotion())
	gsl_odeiv2_driver *odeDriver;	///< A wrapper for the integration method

	/// @}

	StarParam initialValues;		///< Used to read and hold the observational values from file @a StarParam::fromFile

	const OrbitType orbitType;		///< Which orbit calculation method to use

	///< The method will also determine whether #cartesian are the
	///< coordinates in an inertial system or a rotating system.
	///< If \a orbitType is eiher #Numeric or #Linear, the coordinate
	///< system is an inertial system that does not participate
	///< in the galactic rotation. The origin of the system is the
	///< Sun at the epoche #time=0. This means that the Sun is
	///< moving out of the origin with ongoing time.
	///<
	///< If \a orbitType is #Epicycle, the coordinate system is
	///< rotating around the galactic center and the origin is
	///< always the Sun.
	///<
	///< \sa initUVW()

	Vector<mytype> initLB();
	void initXYZ();
	void initUVW();
	void initSigmaXYZ();
	void initSigmaUVW();
	static Matrix<mytype> initT();
	Matrix<mytype> B( mytype alpha, mytype delta ) const;
	static int eqnOfMotion(double, const double y[], double dydt[], void *params);
	static int jacEqnOfMotion(double, const double y[], double *dfdy, double dfdt[], void *params);
	void initODE();
	int numericOrbit(const mytype t);

public:
	StarGalacticParam(OrbitType orbit, bool isSun=false, mytype epoch=0.0);
	StarGalacticParam(QString filename, int linenumber, OrbitType orbit);
	virtual ~StarGalacticParam();
	void initValues(QString filename, int linenumber);
	static Vector<mytype> lbRFromXYZ(const Vector<mytype> xyz, bool allowNegativeLongitude = false);
	static Vector<mytype> equatorialFromGalactic(const Vector<mytype> gal);
	static Vector<mytype> galacticFromEquatorial(const Vector<mytype> equ);
	static Vector<mytype> equPMFromGalPM(Vector<mytype> gal, mytype alpha, mytype delta);
	static Vector<mytype> equPMFromEclPM(Vector<mytype> ecl, mytype alpha, mytype delta);
	Vector<mytype> pmFromVelocities(mytype alpha, mytype delta, mytype para, mytype u, mytype v, mytype w) const;
	mytype varyStartValue(StarParam::Uncertainty what, Rng &rng);
	bool varyAllStartValues(Rng &rng);
	void coordsFromStartValues(bool includeSigmas);

	/// @name Return individual components
	/// @{
	mytype getX() const {return cartesian[0];}				///< [pc]
	mytype getY() const {return cartesian[1];}				///< [pc]
	mytype getZ() const {return cartesian[2];}				///< [pc]
	mytype getU() const {return cartesian[3];}				///< [pc/yr]
	mytype getV() const {return cartesian[4];}				///< [pc/yr]
	mytype getW() const {return cartesian[5];}				///< [pc/yr]
	mytype getRadius() const {return initialValues.getRadius();}	///< [pc]
	QString getHIP() const {return initialValues.getHIP();}
	QString getRem() const {return initialValues.getRem();}
	mytype getTime() const {return time;}
	bool isAssoc() const {return initialValues.isAssoc();}
	unsigned int getLineNumber() const {return initialValues.getLineNumber();}
	QString getInputFileName() const {return initialValues.getInputFileName();}

	/// @}

	/// @name Starting values for the orbit
	/// @{
	mytype getInitialRA() const {return initialValues.getRA();}			///< [°]
	mytype getInitialDec() const {return initialValues.getDec();}		///< [°]
	mytype getInitialPara() const {return initialValues.getParallax();}	///< [mas]
	mytype getInitialU() const {return initialValues.getU(); }	///< [km/s]
	mytype getInitialV() const {return initialValues.getV(); }	///< [km/s]
	mytype getInitialW() const {return initialValues.getW(); }	///< [km/s]

	/// @}

	void zeroSigmas();
	void initStartingValues();
	OrbitType getOrbitType() const {return orbitType;} ///< Return the orbit type used
	int epicyclicOrbit(const mytype t, bool updateErrors, bool useMaxError=false);
	void epicyclicError(const mytype t, bool maxError=false);
	int linearOrbit(const mytype t, bool updateErrors, bool useMaxErrors=false);
	void linearError(const mytype t, bool maxError=false);
	int orbitAt(const mytype t);
	int progressOrbit(const mytype t);
	void resetODE(const double hstart);
	Matrix<mytype> coordsHeliocentric() const;
	Matrix<mytype> coordsInertial() const;
	mytype getDistance(bool heliocentric=true) const;
	mytype getSpaceVelocity(bool heliocentric) const;
	mytype getRadialVelocity(bool heliocentric) const;
	Vector<mytype> getEquatorial(bool todaysFrame=false) const;
	Vector<mytype> getGalactic(bool todaysFrame=false, bool allowNegativeLongitude = false) const;
	static void printConstants(QTextStream &out);

	/// @name Some formatted output of current values
	/// @{
	void printInfo(QTextStream &out) const;
	void printStartValues(QTextStream &out, bool prependHash) const;
	void printEquatorial(QTextStream &out, bool heliocentric=true) const;
	void printGalactic(QTextStream &out, bool heliocentric=true) const;
	void printXYZ(QTextStream &out, bool heliocentric) const;
	void printUVW(QTextStream &out, bool heliocentric) const;
	void printPM(QTextStream &out, bool heliocentric=true) const;
	void printSigmas(QTextStream &out) const;
	/// @}

	mytype operator ||(const StarGalacticParam &st) const;
	mytype distanceTo(const StarGalacticParam &st, char direction) const;
	bool isCloseTo(const StarGalacticParam &st, double withinSigma=1.) const;
	/**
	 * @brief Convert velocities from km/s to pc/yr
	 * @param kms Velocity in km/s
	 * @return Velocity in pc/yr
	 *
	 * 1 pc/yr = 977,792.22 km/s
	 */
	static mytype kms2pcyr(mytype kms) {return kms / 9.7779222e5; }
	/**
	 * @brief Convert velocities from pc/yr to km/s
	 * @param pcyr Velocity in pc/yr
	 * @return Velocity in km/s
	 *
	 * @sa kms2pcyr()
	 */
	static mytype pcyr2kms(mytype pcyr) {return pcyr * 9.7779222e5; }
#define KKK 4.74057
	/*!
	 * \brief Convert velocities from astronomical units per tropical year to km/s
	 * \param auyr Velocity in au/yr
	 * \return Velocity in km/s
	 *
	 * \sa kms2auyr()
	 */
	static mytype auyr2kms(mytype auyr) {return auyr * KKK;}
	/*!
	 * \brief Convert velocities from km/s to astronomical units per tropical year
	 * \param kms Velocity in km/s
	 * \return Velocity in au/yr
	 *
	 * \sa auyr2kms()
	 */
	static mytype kms2auyr(mytype kms) {return kms / KKK;}
#undef KKK
	/*!
	 * \brief writeConfig - for testing purposes only
	 * \param fName Name of the file where the values should go to
	 */
	static void writeConfig(QString fName);
	/*!
	 * \brief Read configuration from file
	 * \param fName File name
	 */
	static void readConfig(QString fName);
	/*!
	 * \brief Convert a SimulParams.Star to the corresponding file name
	 * \param star String 'fName\#lineNr'
	 * \return fName
	 *
	 * \sa SimulParams.Star1
	 */
	static QString simStarToFname( QString star )
	{
		return star.left( star.lastIndexOf( '#' ) );
	}
	/*!
	 * \brief Convert a SimulParams.Star expression to the corresponding line number
	 * \param star  String 'fName\#lineNr'
	 * \return lineNr
	 *
	 * \sa SimulParams.Star1
	 */
	static int simStarToLineNr( QString star )
	{
		return star.split( '#' ).last().toInt();
	}
	/*!
	 * \brief Convert a SimulParams.Star expression to a list of line numbers
	 * \param star String, <b>included in double quotes</b>, e.g. "fName#1,3,18,4..8,5"
	 * \param uniqueAndSorted bool, remove duplicates and sort values in ascending order
	 * \return List of integers with all line numbers: (1, 3, 18, 4, 5, 6, 7, 8, 5), or
	 *         if \a uniqueAndSorted was true (1, 3, 4, 5, 6, 7, 8, 18)
	 *
	 * \sa simStarToFname()
	 */
	static QList<int> simStarToLineList( QString star, bool uniqueAndSorted=false )
	{
		QList<int> myList;
		QStringList rangeList( star.split( '#' ).last().simplified().split( ',' ) );
		int result;
		bool ok;

		foreach ( QString range, rangeList ) {
			result = range.toInt( &ok );
			if ( ok )
				myList << result;
			else {
				int subRangeStart ( range.split( ".." ).first().toInt() );
				int subRangeEnd( range.split( ".." ).last().toInt() );
				for ( int i=subRangeStart; i<=subRangeEnd; i++ )
					myList << i;
			}
		}

		// Remove duplicates an sort
		if ( uniqueAndSorted ) {
#if ( QT_VERSION < QT_VERSION_CHECK(5,14,0) )
			myList = myList.toSet().toList();
			qSort( myList );
#else
			myList = QSet<int>(myList.begin(), myList.end()).values();
			std::sort( myList.begin(), myList.end() );
#endif
		}

		return myList;
	}
};

#endif // STARPARAM_H
