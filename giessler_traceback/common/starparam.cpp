// $Id: starparam.cpp 143 2022-07-28 10:52:27Z ifg $

#include <QSettings>
#include <cmath>
#include "starparam.h"
#include "VecNorm.h"

// TODO: choose equinox from command line or config file
#if 0
#define B1950
#endif

// constants to be used:

// rotation of sun around galactic center:
// Omega = 220 km/s / 8 kpc = 0.0275 km/(s * pc) = 867240.0 km / (yr * pc) = 0.2810240E-07 / yr

// Oort constants A and B
// Oort B = -12.37 ± 0.64 km s−1 kpc−1 = -0.1264097E-07 / yr
//            ==> - 2 B = 0.2528194E-07 / yr
// Oort A = 14.82 ± 0.84 km s−1 kpc−1 = 0.1514464E-07 / yr

// epicyclic frequency kappa = 0.039 km/(s * pc) = 0.3985431E-07 / yr

// vertical osci frequency nu = 0.074 km/(s * pc) = 0.7567279E-07 /yr

// for tests however, as in Hoogerwerf 2001 (365.25 d/yr):
// Omega = 219.8 km/s at 8.5 kpc ==> Omega = 0.2644337E-07
// Oort B = -12.4 = -0.1268030E-07
// Oort A = 13.5 = 0.1380517E-07
// ==> - 2 B = 0.2536061E-07
// kappa = sqrt(-4 Omega B) = 0.3662294E-07
// nu = sqrt(4 pi G rho) = stays the same, as H01 do not specify rho nor nu

static mytype omega0, OortA, OortB, kappa, nu;

static mytype distGalcenterLSR;				// [pc]
static mytype velocityLSROrbitingGalaxy;	// [km/s]

// Present velocity of the Sun with respect to Local Standard of Rest (LSR)
static mytype vLSRx;						// [km/s]
static mytype vLSRy;						// [km/s]
static mytype vLSRz;						// [km/s]

// Present coordinates of the Sun with respect of the LSR according to: Voigt, Abriss der Astronomie
static mytype lSRx;							// [pc]
static mytype lSRy;							// [pc]
static mytype lSRz;							// [pc]

//static const mytype omega0 = StarGalacticParam::kms2pcyr( velocityLSROrbitingGalaxy ) / distGalcenterLSR;
//static const mytype OortA = omega0/2.;
//static const mytype OortB = -omega0/2.;
//static const mytype kappa = sqrt( -4. * omega0 * OortB );

#define ODESTEPPER gsl_odeiv2_step_rkf45

static GalaxyParams galaxyParams = {
    4.52e-15,						// gravitational constant [pc³/(Msun*yr²)]
    1e11,							// mass disk [Msun]
    3.4e10,							// mass sphere [Msun]
    StarGalacticParam::kms2pcyr( 128.0 ),	// [pc/yr]
    6500.0,							// [pc]
    260.0,							// [pc]
    700.0,							// [pc]
    12000.0							// [pc]
};

#ifdef B1950
#warning Equinox B1950 is used!
	// B1950.0
	// <a href="http://adsabs.harvard.edu/abs/1987AJ.....93..864J">
	// (Johnson et al., 1987)</a>
	static const mytype alphaNGP = StarParam::deg2rad( 192.25 );
	static const mytype deltaNGP = StarParam::deg2rad( 27.4 );
	static const mytype theta0 =   StarParam::deg2rad( 123.0 );
#else
	// J2000.0
	// mail von RNE (19.11.2012):
	// hier die aktuellen Konstanten des
	// galaktischen Nordpols fuer J2000.0:

	// Rektaszension alpha = 12h 51m 26.282s
	// Deklination   delta = 27deg 07' 42''

	// und positionswinkel des Himmelsnordpols
	// relativ zum Grosskreis, der den gal. Nordpol
	// mit der gal. Laenge l=0 deg verbindet, ist

	// 122.932 deg

	static const mytype alphaNGP = StarParam::deg2rad( (12. + 51./60. + 26.282/3600.)*15. );
	static const mytype deltaNGP = StarParam::deg2rad( 27. + 7./60. + 42./3600. );
	static const mytype theta0 =   StarParam::deg2rad( 122.932 );
#endif

	// Obliquity of ecliptic (J2000.0)
	// From Hipparcos catalogue, §1.2, table 1.2.2.
	// ε = 23° 26' 21.448"

	static const mytype epsilon = StarParam::deg2rad( 23.4392911111 );
	static const mytype sinepsilon = sin( epsilon );
	static const mytype cosepsilon = cos( epsilon );

SimulParams simulParams;

StarParam::StarParam()
    : RA(0), Dec(0), para(0), para0(0), p_err(0),
      RV(0), RV0(0), RV_err(0), PMa(0), PMa0(0), PMa_err(0),
      PMd(0), PMd0(0), PMd_err(0), radius(0), HIP("Sun"),
      rem(""), distType(Gaussian)
{
}

/**
 * @brief Read parameters from an input file
 * @param file Input file
 * @param lineNumber Line# to read
 * @return Number of bytes read
 */
int StarParam::readValuesFromFile(QFile &file, int lineNumber)
{
	int rc = 0;
	if ( !file.open( QIODevice::ReadOnly ) )
		qFatal( "StarParam::readValuesFromFile(): cannot open file '%s'",
		        file.fileName().toLatin1().data() );

	char line[1024];

	for ( int i=0; i<lineNumber; i++ ) {
		if ( ( rc = file.readLine(line, sizeof(line)) ) < 1 ) {
			qFatal( "StarParam::readValuesFromFile(): no line #%i in '%s'",
			        lineNumber,
			        file.fileName().toLatin1().data() );
		}
	}
	file.close();

	QString qs = QString::fromUtf8( line );
	QStringList paramList = qs.simplified().split( ' ' );

	coordType = static_cast<CoordType>( paramList.at( 0 ).toInt() );
	int minListSize = coordType < 100 ? 0 : 1;
	Vector<mytype> vec(3);

	switch ( coordType ) {
		case Equatorial:
		case EquatorialAssoc:
		// We know that PMa and PMd are parts of unions, so:
		case EquatorialEPM:
		case EquatorialEPMAssoc:
			minListSize += 17;
			if ( paramList.size() < minListSize ) {
				qFatal( "StarParam::readValuesFromFile(): line #%i in '%s' too short",
				        lineNumber,
				        file.fileName().toLatin1().data() );
			}

			vec[0] = paramList.at( 1 ).toDouble();
			vec[1] = paramList.at( 2 ).toDouble();
			vec[2] = paramList.at( 3 ).toDouble();

			RA = hms2deg( vec, paramList.at( 1 ).startsWith( '-' ) );

			vec[0] = paramList.at( 4 ).toDouble();
			vec[1] = paramList.at( 5 ).toDouble();
			vec[2] = paramList.at( 6 ).toDouble();

			Dec = dms2deg( vec, paramList.at( 4 ).startsWith( '-' ) );

			para = para0 = (mytype)   paramList.at(  7 ).toDouble();
			p_err        = (mytype)   paramList.at(  8 ).toDouble();
			RV = RV0     = (mytype)   paramList.at(  9 ).toDouble();
			RV_err       = (mytype)   paramList.at( 10 ).toDouble();
			distType     = (DistType) paramList.at( 11 ).toInt();
			PMa = PMa0   = (mytype)   paramList.at( 12 ).toDouble();
			PMa_err      = (mytype)   paramList.at( 13 ).toDouble();
			PMd = PMd0   = (mytype)   paramList.at( 14 ).toDouble();
			PMd_err      = (mytype)   paramList.at( 15 ).toDouble();
			break;

		case GalacticUVW:
		case GalacticUVWAssoc:
			minListSize += 12;
			if ( paramList.size() < minListSize ) {
				qFatal( "StarParam::readValuesFromFile(): line #%i in '%s' too short",
				        lineNumber,
				        file.fileName().toLatin1().data() );
			}

			l            = (mytype) paramList.at(  1 ).toDouble();
			b            = (mytype) paramList.at(  2 ).toDouble();

			para = para0 = (mytype) paramList.at(  3 ).toDouble();
			p_err        = (mytype) paramList.at(  4 ).toDouble();

			U = U0       = (mytype) paramList.at(  5 ).toDouble();
			U_err        = (mytype) paramList.at(  6 ).toDouble();
			V = V0       = (mytype) paramList.at(  7 ).toDouble();
			V_err        = (mytype) paramList.at(  8 ).toDouble();
			W = W0       = (mytype) paramList.at(  9 ).toDouble();
			W_err        = (mytype) paramList.at( 10 ).toDouble();
			distType      = Gaussian;
			break;

		case GalacticPM:
		case GalacticPMAssoc:
			minListSize += 13;
			if ( paramList.size() < minListSize ) {
				qFatal( "StarParam::readValuesFromFile(): line #%i in '%s' too short",
				        lineNumber,
				        file.fileName().toLatin1().data() );
			}

			l            = (mytype)   paramList.at(  1 ).toDouble();
			b            = (mytype)   paramList.at(  2 ).toDouble();

			para = para0 = (mytype)   paramList.at(  3 ).toDouble();
			p_err        = (mytype)   paramList.at(  4 ).toDouble();
			RV = RV0     = (mytype)   paramList.at(  5 ).toDouble();
			RV_err       = (mytype)   paramList.at(  6 ).toDouble();
			distType     = (DistType) paramList.at(  7 ).toInt();
			PMl = PMl0   = (mytype)   paramList.at(  8 ).toDouble();
			PMl_err      = (mytype)   paramList.at(  9 ).toDouble();
			PMb = PMb0   = (mytype)   paramList.at( 10 ).toDouble();
			PMb_err      = (mytype)   paramList.at( 11 ).toDouble();
			break;

		default:
			qFatal( "%s:#%i: unknown coordType %i\n",
			        file.fileName().toLatin1().data(),
			        lineNumber,
			        coordType );
			break;
	} /*switch( coordType )*/

	/* Associations */
	if ( coordType > 100 )
		radius = paramList.at(minListSize-2).toDouble();

	HIP = paramList.at(minListSize-1);

	if ( paramList.size() > minListSize )
		rem = qs.trimmed().section( QRegExp("\\s+"), minListSize );

	fromFile = file.fileName();
	this->lineNumber = lineNumber;

	/* We know of the union character */
	p_err   = fabs( p_err );
	RV_err  = fabs( RV_err );
	PMa_err = fabs( PMa_err );
	PMd_err = fabs( PMd_err );

	return rc;
}

/*!
 * \brief Vary one of the error-prone attributes within a "reasonable" interval
 * \param what The attribute to vary
 * \param rng Reference to a random number generator
 * \return The new value of the varied attribute
 *
 * This function is used to simulate a measurement of one of the attributes
 * read from the input file. The new value is also saved in the corresponding
 * attribute.
 *
 * In case of neutron stars the radial velocity is usually unknown. It
 * is possible to generate varying velocities by assuming that the space
 * velocities of neutron stars are Maxwell distributed with a 1-D rms
 * @f$\sigma@f$ of 265 km/s.
 *
 * Another possibility is the assumption of a radial velocity uniformly
 * distributed in the intervall -1000 ... +1000 km/s.
 *
 * In all other cases the observables are varied with a Gaussian
 * distribution according to their uncertainties. Parallaxes are limited
 * such that the distances are between 0.1 pc an 25 kpc. The program will barf
 * if it cannot find a parallax with the given parameters. The same holds
 * for radial velocities of a Maxwell distribution.
 *
 * To obtain correct results it is important to vary #RadialVelocity after all
 * other values.
 *
 * \sa
 * <a href="http://adsabs.harvard.edu/abs/2005MNRAS.360..974H">
 * Hobbs et al. (2005): A statistical study of 233 pulsar proper motions, MNRAS 360, 974-992
 * </a>
 *
 * \sa
 * #StarParam::DistType
 */
mytype StarParam::varyValue(StarParam::Uncertainty what, Rng &rng)
{
	const mytype k = StarGalacticParam::auyr2kms( 1. );
	mytype result = 0.0;
	bool offLimits;
	int plxcounter = 0;

	switch ( what ) {
		case Parallax:
			do {
				if ( ++plxcounter > 100 ) {
					qCritical( "varyValue(#%i): failed to find parallax, mean=%g, sigma=%g",
					        getLineNumber(), para0, p_err );
					result = para = nan( "parallax" );
					return result;
				}
				result = para = rng.gaussian( para0, p_err );
				// TODO: how do we limit the parallax???
				// [para] = mas
				// => Limit:            > 0.0001 kpc               < 25 kpc
				offLimits = (result >= (1./.0001)) || (result <= (1./25.));
			} while ( offLimits );
			break;
		case RadialVelocity:
			switch ( distType ) {
				case Gaussian:
					result = RV = rng.gaussian( RV0, RV_err );
					break;
				case Uniform:
					// NOTE: we do not take the transversal velocity into account
					// TODO: lower and upper limit could be given in input file
					result = RV = rng.uniform(-1000.0, +1000.0);
					break;
				case Maxwellian:
					// NOTE: that's how Nina did it
					// TODO: does the Maxwell distribution refer to the corrected
					// or uncorrected space velocities (LSR, sun)?
					mytype vTrans = k/para * sqrt( PMa*PMa + PMd*PMd );
					mytype vSpace;
					int counter = 0;
					int cnt = 0;

					do {
						if ( ++cnt > 100 ) {
							qCritical( "varyValue(#%i): failed to find Maxwellian RV, para=%g, PMa=%g, PMd=%g, v_trans=%g, v_trans0=%g",
							           getLineNumber(),para, PMa, PMd, vTrans, (k/para0 * sqrt(PMa0*PMa0 + PMd0*PMd0)) );
							result = RV = nan( "radial velocity" );
							return result;
						}
						// Try five times and then shuffle again
						if ( ++counter > 5 ) {
							counter = 0;
							varyValue( Parallax, rng );
							varyValue( PMAlpha,  rng );
							varyValue( PMDelta,  rng );
							vTrans = k/para * sqrt( PMa*PMa + PMd*PMd );
						}
						// Nina used the equivalent of 250.66 (vmean = 400)
						vSpace = rng.maxwellian( 265. );
					} while ( vSpace < vTrans );

					// TODO: taking transversal velocity into account will distort Maxwellian distribution
					result = RV = sqrt( vSpace*vSpace - vTrans*vTrans ) * rng.randomsign();
					break;
			}
			break;
		case PMAlpha:
			result = PMa = rng.gaussian( PMa0, PMa_err );
			break;
		case PMDelta:
			result = PMd = rng.gaussian( PMd0, PMd_err );
			break;
		case PMLong:
			result = PMl = rng.gaussian( PMl0, PMl_err );
			break;
		case PMLat:
			result = PMb = rng.gaussian( PMb0, PMb_err );
			break;
		case PMELong:
			result = PMlambda = rng.gaussian( PMlambda0, PMlambda_err );
			break;
		case PMELat:
			result = PMbeta = rng.gaussian( PMbeta0, PMbeta_err );
			break;
		case VelU:
			result = U = rng.gaussian( U0, U_err );
			break;
		case VelV:
			result = V = rng.gaussian( V0, V_err );
			break;
		case VelW:
			result = W = rng.gaussian( W0, W_err );
			break;
	}

	return result;
}

/*!
 * \brief Vary all error-prone attributes (#para, #RV, #PMa, and #PMd)
 * \param rng Reference to a random number generator
 */
bool StarParam::varyAllValues(Rng &rng)
{
	bool rc = true;
	switch ( coordType ) {
		case Equatorial:
		case EquatorialAssoc:
			varyValue( Parallax, rng );
			varyValue( PMAlpha, rng );
			varyValue( PMDelta, rng );
			// This might refer to a neutron star, so RadialVelocity needs to be last
			varyValue( RadialVelocity, rng );
			break;
		case GalacticUVW:
		case GalacticUVWAssoc:
			varyValue( Parallax, rng );
			varyValue( VelU, rng );
			varyValue( VelV, rng );
			varyValue( VelW, rng );
			break;
		case GalacticPM:
		case GalacticPMAssoc:
			varyValue( Parallax, rng );
			varyValue( PMLong, rng );
			varyValue( PMLat, rng );
			// This might refer to a neutron star, so RadialVelocity needs to be last
			varyValue( RadialVelocity, rng );
			break;
		case EquatorialEPM:
		case EquatorialEPMAssoc:
			varyValue( Parallax, rng );
			varyValue( PMELong, rng );
			varyValue( PMELat, rng );
			// This might refer to a neutron star, so RadialVelocity needs to be last
			varyValue( RadialVelocity, rng );
			break;
	} /*switch ( coordType )*/
	if ( isnan(para) || isnan(RV) ) rc = false;
	return rc;
}

#if ( QT_VERSION < QT_VERSION_CHECK(5,14,0) )
# define ENDL endl
#else
# define ENDL Qt::endl
#endif

/**
 * @brief Print the values read from the input file
 * for debugging purposes
 */
void StarParam::printValues(QTextStream &out, bool prependHash) const
{
	const QString hash   = prependHash ? "# " : "";

	out << ENDL;
	out <<         QString("%1CoordType = %2")
	               .arg(hash)
	               .arg(coordType) << ENDL;

	QString radVel = QString( " (but use %1 distribution instead)")
	                 .arg( distType==Maxwellian ? "Maxwell" : "Uniform" );
	switch ( coordType ) {
		case Equatorial:
		case EquatorialAssoc:
			out << QString::fromUtf8( "%1α         = %2")
			       .arg(hash)
			       .arg(deg2hms( RA )) << ENDL;
			out << QString::fromUtf8( "%1δ         = %2")
			       .arg(hash)
			       .arg(deg2dms( Dec )) << ENDL;
			out << QString::fromUtf8( "%1π         = (%L2 ± %L3) mas (%L4, %L5…%L6 pc)" )
			       .arg(hash)
			       .arg( (double) para )
			       .arg( (double) p_err )
			       .arg( 1000./para, 0, 'f', 0 )
			       .arg( 1000./(para+p_err), 0, 'f', 0 )
			       .arg( 1000./(para-p_err), 0, 'f', 0 ) << ENDL;
			out << QString::fromUtf8( "%1RV        = (%L2 ± %L3) km/s%4" )
			       .arg(hash)
			       .arg( (double) RV )
			       .arg( (double) RV_err )
			       .arg( distType==Gaussian ? "" : radVel ) << ENDL;
			out << QString::fromUtf8( "%1µ_α*      = (%L2 ± %L3) mas/yr" )
			       .arg(hash)
			       .arg( (double) PMa )
			       .arg( (double) PMa_err ) << ENDL;
			out << QString::fromUtf8( "%1µ_δ       = (%L2 ± %L3) mas/yr" )
			       .arg(hash)
			       .arg( (double) PMd )
			       .arg( (double) PMd_err ) << ENDL;
			break;
		case GalacticUVW:
		case GalacticUVWAssoc:
			out << QString::fromUtf8( "%1l         = %L2°")
			       .arg(hash)
			       .arg((double) l, 0, 'f') << ENDL;
			out << QString::fromUtf8( "%1b         = %L2°")
			       .arg(hash)
			       .arg((double) b, 0, 'f') << ENDL;
			out << QString::fromUtf8( "%1π         = (%L2 ± %L3) mas (%L4, %L5…%L6 pc)" )
			       .arg(hash)
			       .arg( (double) para )
			       .arg( (double) p_err )
			       .arg( 1000./para, 0, 'f', 0 )
			       .arg( 1000./(para+p_err), 0, 'f', 0 )
			       .arg( 1000./(para-p_err), 0, 'f', 0 ) << ENDL;
			out << QString::fromUtf8( "%1U         = (%L2 ± %L3) km/s" )
			       .arg(hash)
			       .arg( (double) U )
			       .arg( (double) U_err ) << ENDL;
			out << QString::fromUtf8( "%1V         = (%L2 ± %L3) km/s" )
			       .arg(hash)
			       .arg( (double) V )
			       .arg( (double) V_err ) << ENDL;
			out << QString::fromUtf8( "%1W         = (%L2 ± %L3) km/s" )
			       .arg(hash)
			       .arg( (double) W )
			       .arg( (double) W_err ) << ENDL;
			break;
		case GalacticPM:
		case GalacticPMAssoc:
			out << QString::fromUtf8( "%1l         = %L2°")
			       .arg(hash)
			       .arg((double) l, 0, 'f') << ENDL;
			out << QString::fromUtf8( "%1b         = %L2°")
			       .arg(hash)
			       .arg((double) b, 0, 'f') << ENDL;
			out << QString::fromUtf8( "%1π         = (%L2 ± %L3) mas (%L4, %L5…%L6 pc)" )
			       .arg(hash)
			       .arg( (double) para )
			       .arg( (double) p_err )
			       .arg( 1000./para, 0, 'f', 0 )
			       .arg( 1000./(para+p_err), 0, 'f', 0 )
			       .arg( 1000./(para-p_err), 0, 'f', 0 ) << ENDL;
			out << QString::fromUtf8( "%1RV        = (%L2 ± %L3) km/s%4" )
			       .arg(hash)
			       .arg( (double) RV )
			       .arg( (double) RV_err )
			       .arg( distType==Gaussian ? "" : radVel ) << ENDL;
			out << QString::fromUtf8( "%1µ_l*      = (%L2 ± %L3) mas/yr" )
			       .arg(hash)
			       .arg( (double) PMl )
			       .arg( (double) PMl_err ) << ENDL;
			out << QString::fromUtf8( "%1µ_b       = (%L2 ± %L3) mas/yr" )
			       .arg(hash)
			       .arg( (double) PMb )
			       .arg( (double) PMb_err ) << ENDL;
			break;
		case EquatorialEPM:
		case EquatorialEPMAssoc:
			out << QString::fromUtf8( "%1α         = %2")
			       .arg(hash)
			       .arg(deg2hms( RA )) << ENDL;
			out << QString::fromUtf8( "%1δ         = %2")
			       .arg(hash)
			       .arg(deg2dms( Dec )) << ENDL;
			out << QString::fromUtf8( "%1π         = (%L2 ± %L3) mas (%L4, %L5…%L6 pc)" )
			       .arg(hash)
			       .arg( (double) para )
			       .arg( (double) p_err )
			       .arg( 1000./para, 0, 'f', 0 )
			       .arg( 1000./(para+p_err), 0, 'f', 0 )
			       .arg( 1000./(para-p_err), 0, 'f', 0 ) << ENDL;
			out << QString::fromUtf8( "%1RV        = (%L2 ± %L3) km/s%4" )
			       .arg(hash)
			       .arg( (double) RV )
			       .arg( (double) RV_err )
			       .arg( distType==Gaussian ? "" : radVel ) << ENDL;
			out << QString::fromUtf8( "%1µ_λ*      = (%L2 ± %L3) mas/yr" )
			       .arg(hash)
			       .arg( (double) PMlambda )
			       .arg( (double) PMlambda_err ) << ENDL;
			out << QString::fromUtf8( "%1µ_β       = (%L2 ± %L3) mas/yr" )
			       .arg(hash)
			       .arg( (double) PMbeta )
			       .arg( (double) PMbeta_err ) << ENDL;
			break;
	} /*switch ( valueType )*/

	/* Associations */
	if ( coordType > 100 ) {
		out <<     QString("%1Radius    = %2 pc")
		           .arg(hash)
		           .arg( radius ) << ENDL;
	}
	out <<         QString("%1ID        = %2")
	               .arg(hash)
	               .arg( HIP ) << ENDL;
	out <<         QString("%1Comment   = %2")
	               .arg(hash)
	               .arg( rem ) << ENDL;
}

/**
 * @brief Convert angles from degrees to radian
 * @param degrees Angle in degrees to convert
 * @return Angle in radian
 */
mytype StarParam::deg2rad(mytype degrees)
{
	return degrees / 180. * M_PI;
}

/**
 * @brief Convert angles from radian to degrees
 * @param radian Angle in radian to convert
 * @return Angle in degrees
 */
mytype StarParam::rad2deg(mytype radian)
{
	return radian / M_PI * 180.;
}

/**
 * @brief Convert angle from decimal degrees to deg, arcmin, arcsec
 * @param degrees Floating point angle to convert [°]
 * @return String representation of @a degrees
 */
QString StarParam::deg2dms(mytype degrees)
{
	Vector<mytype> dms(3);

	mytype fraction = modf( fabs( degrees ), &dms[0] );
	dms[2] = modf( fraction * 60.0, &dms[1] ) * 60.0;

	/* Catch rounding errors */
	if ( QString("%1").arg( dms[2] ) == "60" ) {
		dms[1] += 1.0;
		dms[2] = 0.0;
	}

	/* dto. */
	if ( QString("%1").arg( dms[1] ) == "60" ) {
		dms[0] += 1.0;
		dms[1] = 0.0;
	}

	return QString( "%1%2%3 %4%5 %L6%7" )
	        .arg( degrees<0 ? "-" : "" )
	        .arg( dms[0] )
	        .arg( QString::fromUtf8("°") )
	        .arg( dms[1] )
	        .arg( QString::fromUtf8("′") )
	        .arg( dms[2], 0, 'f' )
	        .arg( QString::fromUtf8( "″") );
}

/**
 * @brief Convert angle from decimal degrees to hour, minutes, seconds
 * @param degrees Floating point angle to convert [°]
 * @return String representation of @a degrees
 */
QString StarParam::deg2hms(mytype degrees)
{
	Vector<mytype> hms(3);

	mytype fraction = modf( fabs( degrees/15.0 ), &hms[0] );
	hms[2] = modf( fraction * 60.0, &hms[1] ) * 60.0;

	/* Catch rounding errors */
	if ( QString("%1").arg( hms[2] ) == "60" ) {
		hms[1] += 1.0;
		hms[2] = 0.0;
	}

	/* dto. */
	if ( QString("%1").arg( hms[1] ) == "60" ) {
		hms[0] += 1.0;
		hms[1] = 0.0;
	}

	return QString( "%1%2h %3m %L4s" )
	        .arg( degrees<0 ? "-" : "" )
	        .arg( hms[0] )
	        .arg( hms[1] )
	        .arg( hms[2], 0, 'f' );
}

/**
 * @brief Convert angle from degrees, arcmin, arcsec to decimal degrees
 * @param dms %Vector containing degrees, arcmin, and arcsec
 * @param isNegative Whether the angle is smaller than 0
 * @return Decimal value [°]
 */
mytype StarParam::dms2deg(Vector<mytype> dms, bool isNegative)
{
	const mytype deg = fabs( dms[0] ) + (dms[1]/60.0) + (dms[2]/3600.0);

	return isNegative ? -deg : deg;
}

/**
 * @brief Convert angle from hours, minutes, seconds to decimal degrees
 * @param hms %Vector containing hours, minutes, and seconds
 * @param isNegative Whether the angle is smaller than 0
 * @return Decimal value [°]
 */
mytype StarParam::hms2deg(Vector<mytype> hms, bool isNegative)
{
	const mytype deg = dms2deg( hms, isNegative );

	return deg * 15.0;
}

/**
 * @brief Constructor
 * @param orbit Calculation method for the orbit
 * @param isSun Whether this will be representing the Sun
 * @param epoch Calculate orbit for given time [yr] (for Sun only)
 *
 * Set all attributes to default values (zero) or create an
 * instance representing our Sun, and
 * initialize transformation matrix @b T and constant @a k
 */
StarGalacticParam::StarGalacticParam(OrbitType orbit, bool isSun, mytype epoch)
    :time(0.0),
      x0(0),  y0(0), z0(0), u0(0), v0(0), w0(0),
      sigmaX0(0), sigmaY0(0), sigmaZ0(0),
      sigmaX(0), sigmaY(0), sigmaZ(0),
      sigmaU0(0), sigmaV0(0), sigmaW0(0),
      sigmaU(0), sigmaV(0), sigmaW(0),
      T(initT()), orbitType(orbit)
{
	if ( isSun ) {
		if ( orbitType != Epicycle ) {
			x0 = x() = 0.0;
			y0 = y() = 0.0;
			z0 = z() = 0.0;
			u0 = u() = kms2pcyr( vLSRx );
			v0 = v() = kms2pcyr( vLSRy + velocityLSROrbitingGalaxy );
			w0 = w() = kms2pcyr( vLSRz );
		} else {
			/* for epicyclic orbits */
			x0 = x() = 0.0;
			y0 = y() = 0.0;
			z0 = z() = 0.0;
			u0 = u() = /*kms2pcyr( vLSRx )*/0.;
			v0 = v() = /*kms2pcyr( vLSRy )*/0.;
			w0 = w() = /*kms2pcyr( vLSRz )*/0.;
		}
	} else /* it is not Sun */
		for ( int i=0; i<6; i++) cartesian[i] = 0.0;

	if ( orbitType == Numeric )
		initODE();

	if ( epoch != 0.0 ) {
		assert( isSun );
		if ( (epoch > 0.0 && simulParams.StepSize < 0.0) ||
		     (epoch < 0.0 && simulParams.StepSize > 0.0) )
			resetODE( -simulParams.StepSize );

		orbitAt( epoch );
	}
}

/**
 * @brief Constructor
 * @param filename File with input parameters
 * @param linenumber Which line in @a filename to read
 * @param orbit Calculation method for the orbit
 */
StarGalacticParam::StarGalacticParam(QString filename, int linenumber, OrbitType orbit)
    : time(0.0), T(initT()), orbitType(orbit)
{
	initValues( filename, linenumber );
	if ( orbitType == Numeric )
		initODE();
}

StarGalacticParam::~StarGalacticParam()
{
	if ( orbitType == Numeric )
		gsl_odeiv2_driver_free( odeDriver );
}

/**
 * @brief Initialize members
 * @param filename File with input parameters
 * @param linenumber Which line in @a filename to read
 */
void StarGalacticParam::initValues(QString filename, int linenumber)
{
	QFile qFile( filename );

	initialValues.readValuesFromFile( qFile, linenumber );

	coordsFromStartValues( true );
	time = 0.;
}

/**
 * @brief Calculate galactic coordinates from cartesian ones
 * @param xyz Vector with \b heliocentric(!) cartesian components [pc]
 * @param allowNegativeLongitude If \a true, Gal. long. will be in the range -180°...+180°,
 *        otherwise 0°...360°
 * @return Vector with l, b, and R in degrees and pc, resp.
 */
Vector<mytype> StarGalacticParam::lbRFromXYZ(const Vector<mytype> xyz, bool allowNegativeLongitude)
{
	const mytype
	        x = xyz[0],
	        y = xyz[1],
	        z = xyz[2];

	const mytype r = sqrt( x*x + y*y + z*z );
	const mytype theta = acos( z/r );
	const mytype phi = atan2( y, x );

	Vector<mytype> gal(3);

	if ( allowNegativeLongitude ) {
		gal[0] = StarParam::rad2deg( phi );
	} else {
		gal[0] = StarParam::rad2deg( phi > 0 ? phi : 2.*M_PI + phi );	// gal. longitude l
	}
	gal[1] = 90. - StarParam::rad2deg( theta );							// gal. latitude b
	gal[2] = r;															// distance

	return gal;
}

/**
 * @brief Calculate equatorial coordinates from galactic ones
 * @param gal Vector with l, b, and R in degrees and pc, resp.
 * @return Vector with alpha, delta, and R in degrees and pc, resp.
 *
 * \sa
 * <a href="http://adsabs.harvard.edu/abs/1987AJ.....93..864J">
 * Johnson et al., 1987</a>
 */
Vector<mytype> StarGalacticParam::equatorialFromGalactic(const Vector<mytype> gal)
{
	const mytype l = StarParam::deg2rad( gal[0] );
	const mytype b = StarParam::deg2rad( gal[1] );

	Vector<mytype> equ(3), gal1(3), result(3);

	gal1[0] = cos( b ) * cos( l );
	gal1[1] = cos( b ) * sin( l );
	gal1[2] = sin( b );

	// TODO: check whether T is always a rotation (orthogonal) matrix
	// in which case we could use the much faster transpose() instead
	// of inverse()
	equ = initT().inverse() * gal1;

	const mytype delta = asin( equ[2] );

	const mytype cosAlpha = equ[0]/* / cos( delta )*/;
	const mytype sinAlpha = equ[1]/* / cos( delta )*/;

	const mytype alpha = atan2( sinAlpha, cosAlpha );

	result[0] = StarParam::rad2deg( alpha >= 0. ? alpha : 2.*M_PI + alpha );
	result[1] = StarParam::rad2deg( delta );
	result[2] = gal[2];

	return result;
}

/*!
 * \brief Calculate galactic coordinates from equatorial ones
 * \param equ Vector with alpha, delta, and R in degrees and pc, resp.
 * \return Vector with l, b, and R in degrees and pc, resp.
 *
 * \sa
 * <a href="http://adsabs.harvard.edu/abs/1987AJ.....93..864J">
 * Johnson et al., 1987</a>
 */
Vector<mytype> StarGalacticParam::galacticFromEquatorial(const Vector<mytype> equ)
{
	const mytype alpha = StarParam::deg2rad( equ[0] );
	const mytype delta = StarParam::deg2rad( equ[1] );

	Vector<mytype> equ1(3), gal(3), result(3);

	equ1[0] = cos( delta ) * cos( alpha );
	equ1[1] = cos( delta ) * sin( alpha );
	equ1[2] = sin( delta );

	gal = initT() * equ1;

	// catch rounding errors
	if ( gal[2] > 1. ) gal[2] = 1.;
	else if ( gal[2] < -1. ) gal[2] = -1.;

	const mytype b = asin( gal[2] );

	const mytype cosL = gal[0];
	const mytype sinL = gal[1];

	const mytype l = atan2( sinL, cosL );

	result[0] = StarParam::rad2deg( l > 0. ? l : 2.*M_PI + l );
	result[1] = StarParam::rad2deg( b );
	result[2] = equ[2];

	return result;
}

/**
 * @brief Calculate proper motion in equatorial coordinates from galactic ones
 * @param gal vector containing
 * @f$
 * \left( \!\! \begin{array}{c}
 * \mu_{\ell}^{\star} \\ \mu_{b}
 * \end{array} \!\!\! \right)
 * @f$
 *        in convenient units
 * @param alpha Right ascension of object in rad
 * @param delta Declination of object in rad
 * @return vector containing
 * @f$
 * \left( \!\! \begin{array}{c}
 * \mu_{\alpha}^{\star} \\ \mu_{\delta}
 * \end{array} \!\!\! \right)
 * @f$
 *         in the same units as \a gal
 *
 * The proper motions in galactic longitude and right ascension are corrected
 * for latitude and declination, respectively i.e.,
 * @f$\mu_{\ell}^{\star} = \mu_\ell \cos b@f$
 * and
 * @f$\mu_{\alpha}^{\star} = \mu_{\alpha} \cos \delta @f$
 *
 * The function is not really used in the program but might come in handy from time to time.
 *
 * \sa <a href="http://adsabs.harvard.edu/abs/2013arXiv1306.2945P">http://adsabs.harvard.edu/abs/2013arXiv1306.2945P</a>
 * \sa StarGalacticParam::equPMFromEclPM()
 */
Vector<mytype> StarGalacticParam::equPMFromGalPM(Vector<mytype> gal, mytype alpha, mytype delta)
{
	Vector<mytype> result(2);
	Matrix<mytype> C(2,2);

	const mytype C1 = sin(deltaNGP)*cos(delta) - cos(deltaNGP)*sin(delta)*cos(alpha-alphaNGP);
	const mytype C2 = cos(deltaNGP)*sin(alpha-alphaNGP);

	const mytype cosb = sqrt(C1*C1 + C2*C2);

	C[0][0] =  C1; C[0][1] = C2;
	C[1][0] = -C2; C[1][1] = C1;

	C = 1/cosb * C;

	// Since C is a rotation matrix and therefore orthogonal, we use
	// transpose() instead of inverse() because it's faster
	result = C.transpose() * gal;

	return result;
}

/**
 * @brief Calculate proper motion in equatorial coordinates from ecliptic ones
 * @param ecl vector containing
 * @f$
 * \left( \!\! \begin{array}{c}
 * \mu_{\lambda}^{\star} \\ \mu_{\beta}
 * \end{array} \!\!\! \right)
 * @f$
 *        in convenient units
 * @param alpha Right ascension of object in rad
 * @param delta Declination of object in rad
 * @return vector containing
 * @f$
 * \left( \!\! \begin{array}{c}
 * \mu_{\alpha}^{\star} \\ \mu_{\delta}
 * \end{array} \!\!\! \right)
 * @f$
 *         in the same units as \a ecl
 *
 * The proper motions in ecliptic longitude and right ascension are corrected
 * for latitude and declination, respectively i.e.,
 * @f$\mu_{\lambda}^{\star} = \mu_\lambda \cos \beta@f$
 * and
 * @f$\mu_{\alpha}^{\star} = \mu_{\alpha} \cos \delta @f$
 *
 * The rotation matrix @f$\bf{E}@f$ used in the transformation
 *
 * @f$
 * \left( \!\! \begin{array}{c}
 * \mu_{\lambda}^{\star} \\ \mu_{\beta}
 * \end{array} \!\!\! \right)
 * = \bf{E}
 * \left( \!\! \begin{array}{c}
 * \mu_{\alpha}^{\star} \\ \mu_{\delta}
 * \end{array} \!\!\! \right)
 * @f$
 *
 * is obtained by calculating the time derivative of
 *
 * @f$
 * \displaystyle
 * \sin \beta = \sin \delta \cos \epsilon - \cos \delta \sin \epsilon \sin \alpha
 * ~\left| \,\frac{\textrm{d}}{\textrm{d} t} \right.
 * @f$
 *
 * keeping in mind that
 *
 * @f$
 * \dot{\beta} \equiv \mu_\beta,~ \dot{\alpha} \equiv \mu_\alpha, \dot{\delta} \equiv \mu_\delta
 * @f$
 *
 * This will give the elements
 *
 * @f$
 * \begin{array}{rcl}
 * \cos \beta \cdot E_{21} & = & \sin \epsilon \cos \alpha \\
 * \cos \beta \cdot E_{22} & = & \cos \epsilon \cos \delta + \sin \epsilon \sin \delta \sin \alpha
 * \end{array}
 * @f$
 *
 * The procedure how to get the other elements of the
 * matrix is described in the arXive paper below.
 *
 * \sa <a href="http://adsabs.harvard.edu/abs/2013arXiv1306.2945P">http://adsabs.harvard.edu/abs/2013arXiv1306.2945P</a>
 * \sa StarGalacticParam::equPMFromGalPM()
 */
Vector<mytype> StarGalacticParam::equPMFromEclPM(Vector<mytype> ecl, mytype alpha, mytype delta)
{
	Vector<mytype> result(2);
	Matrix<mytype> F(2,2);

	const mytype F1 = cosepsilon*cos(delta) + sinepsilon*sin(delta)*sin(alpha);
	const mytype F2 = sinepsilon*cos(alpha);

	const mytype cosbeta = sqrt(F1*F1 + F2*F2);

	F[0][0] =  F1; F[0][1] = F2;
	F[1][0] = -F2; F[1][1] = F1;

	F = 1/cosbeta * F;

	// Since F is a rotation matrix and therefore orthogonal, we use
	// transpose() instead of inverse() because it's faster
	result = F.transpose() * ecl;

	return result;
}

/**
 * @brief Calculate radial velocity and proper motion from space velocities
 * @param alpha Right ascension [rad]
 * @param delta Declination [rad]
 * @param para Parallax [arcsec]
 * @param u X component of heliocentric space velocity [km/s]
 * @param v Y component of heliocentric space velocity [km/s]
 * @param w Z component of heliocentric space velocity [km/s]
 * @return Vector with @f$ \rho, \mu^*_{\alpha}, \mu_{\delta}@f$, in km/s and arcsec/yr, resp.
 *
 * The formula that is used for the calculation is mentioned in StarGalacticParam::B().
 *
 * When calling this function, make sure \a u, \a v, and \a w are heliocentric!
 */
Vector<mytype> StarGalacticParam::pmFromVelocities(mytype alpha, mytype delta, mytype para, mytype u, mytype v, mytype w) const
{
	const mytype k = auyr2kms( 1. );
	Vector<mytype> vel(3), pm(3), result(3);

	vel[0] = u;
	vel[1] = v;
	vel[2] = w;

	// TODO: check whether B is always a rotation (orthogonal) matrix
	// in which case we could use the much faster transpose() instead
	// of inverse()
	pm = B(alpha, delta).inverse() * vel;

	result[0] = pm[0];				// radial velocity in km/s
	result[1] = pm[1] * para / k;	// proper motion in RA, corrected for declination, in arcsec/yr
	result[2] = pm[2] * para / k;	// proper motion in declination, in arcsec/yr

	return result;
}

/**
 * @brief Vary one of the error-prone attributes in #initialValues
 * @param what The attribute to vary
 * @param rng Reference to a random number generator
 * @return The new value of the attribute
 *
 * The function will modify only #initialValues. To update the galactic
 * coordinates #x, #y, #z, #u, #v, #w one needs to call
 * coordsFromStartValues() which is time consuming and therefore not
 * called automatically. The returned value might be used e.g., to
 * fill a table prior to calculations.
 *
 * In case of neutron stars make sure that #StarParam::PMAlpha and
 * #StarParam::PMDelta are changed \b before #StarParam::RadialVelocity.
 */
mytype StarGalacticParam::varyStartValue(StarParam::Uncertainty what, Rng &rng)
{
	return initialValues.varyValue(what, rng);
}

/**
 * @brief Vary all error-prone attributes in #initialValues and reset #time
 * @param rng Reference to a random number generator
 * @return @a true if successful
 *
 * The function will modify only #initialValues and #time. To update the galactic
 * coordinates #x, #y, #z, #u, #v, #w one needs to call
 * coordsFromStartValues() which is time consuming and therefore not
 * called automatically.
 */
bool StarGalacticParam::varyAllStartValues(Rng &rng)
{
	time = 0.0;
	return initialValues.varyAllValues( rng );
}

/*!
 * \brief Update attributes after #initialValues has changed
 * \param includeSigmas Whether to calculate also the uncertainties
 *
 * This will calculate the galactic and cartesian coordinates
 * from #initialValues
 */
void StarGalacticParam::coordsFromStartValues(bool includeSigmas)
{
	initXYZ();
	initUVW();
	if ( includeSigmas ) {
		initSigmaXYZ();
		initSigmaUVW();
	}
	initStartingValues();
}

/*!
 * \brief Calculate the distance from the origin of the coordinate system
 * \param heliocentric If \b true, distance to Sun will be returned
 * \return Distance in parsec
 */
mytype StarGalacticParam::getDistance(bool heliocentric) const
{
	const Vector<mytype> coords(heliocentric ? coordsHeliocentric()[0] : coordsInertial()[0] );

	return norm2( coords );
}

/*!
 * \brief Calulate space velocity
 * \param heliocentric Relative to Sun
 * \return Space velocity in km/s
 */
mytype StarGalacticParam::getSpaceVelocity(bool heliocentric) const
{
	const Vector<mytype> coords(heliocentric ? coordsHeliocentric()[1] : coordsInertial()[1] );

	return pcyr2kms( norm2( coords ) );
}

/*!
 * \brief Calculate radial velocity of the object
 * \param heliocentric If \a true RV is always relative to the Sun
 * \return Radial component of space velocity [km/s]
 */
mytype StarGalacticParam::getRadialVelocity(bool heliocentric) const
{
	const Matrix<mytype> coords(heliocentric ? coordsHeliocentric() : coordsInertial());

	return pcyr2kms( normalize(coords[0]) | coords[1] );
}

/*!
 * \brief Get equatorial coordinates
 * \param todaysFrame If \a true, coordinates are relative to the inertial system at time 0.
 *        This is adventurous.
 * \return Vector with components RA [°], Dec [°], and distance [pc]
 *
 * Be aware that earth's precession is \b not taken into account when
 * \a todaysFrame is false. It might be better to use getGalactic() instead.
 */
Vector<mytype> StarGalacticParam::getEquatorial(bool todaysFrame) const
{
	const Vector<mytype> xyz( todaysFrame ? coordsInertial()[0] : coordsHeliocentric()[0] );

	return equatorialFromGalactic( lbRFromXYZ( xyz ) );
}

/*!
 * \brief Get galactic coordinates
 * \param todaysFrame If \a true, coordinates are relative to the inertial system at time 0.
 *        Some people call it 'present day's coordinates'. This is adventurous.
 * \param allowNegativeLongitude If \a true, Galactic longitude will be in the range -180°...+180°,
 *        otherwise 0°...+360°
 * \return Vector with components galactic longitude [°], latitude [°], and distance [pc].
 */
Vector<mytype> StarGalacticParam::getGalactic(bool todaysFrame, bool allowNegativeLongitude) const
{
	const Vector<mytype> xyz(todaysFrame ? coordsInertial()[0] : coordsHeliocentric()[0]);

	return lbRFromXYZ( xyz, allowNegativeLongitude );
}

/**
 * @brief Set uncertainties of current values to zero
 *
 * Usefull for numeric integration
 */
void StarGalacticParam::zeroSigmas()
{
	sigmaX = sigmaY = sigmaZ = 0.;
	sigmaU = sigmaV = sigmaW = 0.;
}

/**
 * @brief Set starting values (x0, y0, ...) to current values
 */
void StarGalacticParam::initStartingValues()
{
	x0 = x(); y0 = y(); z0 = z();
	u0 = u(); v0 = v(); w0 = w();
	sigmaX0 = sigmaX; sigmaY0 = sigmaY; sigmaZ0 = sigmaZ;
	sigmaU0 = sigmaU; sigmaV0 = sigmaV; sigmaW0 = sigmaW;
}

/**
 * @brief Calculate stellar orbit, optionally including errors
 * @param t Time for which to calculate the coordinates, in years
 * @param updateErrors If true, epicyclicError() is called
 * @param useMaxError If true, use maximum error instead of Gaussian error propagation
 *
 * The equations for the epicyclic orbit are given in
 * <a href="http://adsabs.harvard.edu/abs/2006MNRAS.373..993F">
 * Fuchs et al., 2006</a>:
 *
 * @f$
 * \renewcommand{\arraystretch}{2.5}
 * \begin{array}{ccl}
 * X(t) & = & \displaystyle X(0) - \frac{V(0)}{-2B} [1 - \cos(\kappa \,t)]
 *                          + \frac{U(0)}{\kappa} \sin(\kappa \,t), \\
 *
 * U(t) & = & \displaystyle U(0) \cos(\kappa \,t) - \frac{\kappa}{-2B} V(0) \sin(\kappa \,t), \\
 *
 * Y(t) & = & \displaystyle Y(0) + 2A \left[ X(0) - \frac{V(0)}{-2B} \right] t \\
 * & & + \displaystyle \frac{\Omega_0}{-B \,\kappa} V(0) \sin(\kappa \,t)
 *     + \frac{2\,\Omega_0}{\kappa^2}\, U(0) [1 - \cos(\kappa \,t)],\\
 *
 * V(t) & = & \displaystyle \frac{-2B}{\kappa} U(0) \sin(\kappa \,t) + V(0) \cos(\kappa \,t), \\
 *
 * Z(t) & = & \displaystyle \frac{W(0)}{\nu} \sin(\nu \,t) + Z(0) \cos(\nu \,t), \\
 *
 * W(t) & = & W(0) \cos(\nu \,t) - Z(0) \nu \sin(\nu \,t)
 * \end{array}
 * @f$
 */
int StarGalacticParam::epicyclicOrbit(const mytype t, bool updateErrors, bool useMaxError)
{
	const mytype sinnut = sin( nu * t );
	const mytype cosnut = cos( nu * t );
	const mytype sinkappat = sin( kappa * t );
	const mytype coskappat = cos( kappa * t );

	const mytype A = OortA;
	const mytype B = OortB;

	x() = x0 - v0/(-2.*B) * (1. - coskappat) + u0/kappa * sinkappat;
	u() = u0 * coskappat - kappa/(-2.*B) * v0 * sinkappat;

	y() = y0 + 2.*A * ( x0 - v0/(-2.*B) ) * t
	      + omega0/(-B*kappa) * v0 * sinkappat
	      + 2.*omega0/(kappa*kappa) * u0 * (1. - coskappat);
	v() = (-2.*B)/kappa * u0 * sinkappat + v0 * coskappat;

	z() = w0/nu * sinnut + z0 * cosnut;
	w() = w0 * cosnut - z0 * nu * sinnut;

	if ( updateErrors )
		epicyclicError( t, useMaxError );

	time = t;

	return GSL_SUCCESS;
}

/**
 * @brief Calculate errors in position and velocities
 * @param t Time in years
 * @param maxError If true, use maximum error instead of Gaussian error propagation
 *
 * Will calculate uncertainties @f$ \sigma @f$ of the epicyclic orbit according
 * to Gaussian error propagation
 *
 * @f$ \displaystyle
 * \sigma^2 = \sum_i{ \left( \frac{\partial F}{\partial x_i} \right)^{\!\! 2} \sigma_i^2}
 * @f$ ,
 *
 * or maximum error
 *
 * @f$ \displaystyle
 * \Delta F = \sum_i{ \left| \frac{\partial F}{\partial x_i} \right| \Delta x_i}
 * @f$
 *
 * @sa StarGalacticParam::linearError()
 */
void StarGalacticParam::epicyclicError(const mytype t, bool maxError)
{
	//FIXME: this whole error progession stuff needs to be checked
	const mytype coskappat = cos( kappa * t );
	const mytype sinkappat = sin( kappa * t );
	const mytype sinnut = sin( nu * t );
	const mytype cosnut = cos( nu * t );
	mytype sigmaSquare;

	// Error in x
	mytype dX_dX0 = 1.;
	mytype dX_dV0 = (1. - coskappat ) / (2.*OortB);
	mytype dX_dU0 = sinkappat / kappa;

	if ( maxError ) {
		sigmaX = fabs( dX_dX0 ) * sigmaX0
			   + fabs( dX_dV0 ) * sigmaV0
			   + fabs( dX_dU0 ) * sigmaU0;
	} else {
		sigmaSquare = dX_dX0*dX_dX0 * sigmaX0*sigmaX0
					+ dX_dV0*dX_dV0 * sigmaV0*sigmaV0
					+ dX_dU0*dX_dU0 * sigmaU0*sigmaU0;
		sigmaX = sqrt( sigmaSquare );
	}

	// Error in u
	mytype dU_dU0 = coskappat;
	mytype dU_dV0 = kappa*sinkappat/(2.*OortB);

	if ( maxError ) {
		sigmaU = fabs( dU_dU0 ) * sigmaU0
			   + fabs( dU_dV0 ) * sigmaV0;

	} else {
		sigmaSquare = dU_dU0*dU_dU0 * sigmaU0*sigmaU0
					+ dU_dV0*dU_dV0 * sigmaV0*sigmaV0;
		sigmaU = sqrt( sigmaSquare );
	}

	// Error in y
	mytype dY_dY0 = 1.;
	mytype dY_dX0 = 2. * OortA * t;
	mytype dY_dV0 = OortA*t/OortB - omega0*sinkappat / (OortB*kappa);
	mytype dY_dU0 = (2*omega0*(1 - coskappat)) / (kappa * kappa);

	if ( maxError ) {
		sigmaY = fabs( dY_dY0 ) * sigmaY0
			   + fabs( dY_dX0 ) * sigmaX0
			   + fabs( dY_dV0 ) * sigmaV0
			   + fabs( dY_dU0 ) * sigmaU0;
	} else {
		sigmaSquare = dY_dY0*dY_dY0 * sigmaY0*sigmaY0
					+ dY_dX0*dY_dX0 * sigmaX0*sigmaX0
					+ dY_dV0*dY_dV0 * sigmaV0*sigmaV0
					+ dY_dU0*dY_dU0 * sigmaU0*sigmaU0;
		sigmaY = sqrt( sigmaSquare );
	}

	// Error in v
	mytype dV_dU0 = (-2.*OortB*sinkappat)/kappa;
	mytype dV_dV0 = coskappat;

	if ( maxError ) {
		sigmaV = fabs( dV_dU0 ) * sigmaU0
			   + fabs( dV_dV0 ) * sigmaV0;
	} else {
		sigmaSquare = dV_dU0*dV_dU0 * sigmaU0*sigmaU0
					+ dV_dV0*dV_dV0 * sigmaV0*sigmaV0;
		sigmaV = sqrt( sigmaSquare );
	}

	// Error in z
	mytype dZ_dW0 = sinnut / nu;
	mytype dZ_dZ0 = cosnut;

	if ( maxError ) {
		sigmaZ = fabs( dZ_dW0 ) * sigmaW0
			   + fabs( dZ_dZ0 ) * sigmaZ0;
	} else {
		sigmaSquare = dZ_dW0*dZ_dW0 * sigmaW0*sigmaW0
					+ dZ_dZ0*dZ_dZ0 * sigmaZ0*sigmaZ0;
		sigmaZ = sqrt( sigmaSquare );
	}

	// Error in w
	mytype dW_dW0 = cosnut;
	mytype dW_dZ0 = -nu*sinnut;

	if ( maxError ) {
		sigmaW = fabs( dW_dW0 ) * sigmaW0
			   + fabs( dW_dZ0 ) * sigmaZ0;
	} else {
		sigmaSquare = dW_dW0*dW_dW0 * sigmaW0*sigmaW0
					+ dW_dZ0*dW_dZ0 * sigmaZ0*sigmaZ0;
		sigmaW = sqrt( sigmaSquare );
	}
}

/**
 * @brief Calculate new position of the star assuming a linear motion
 *        determined by its space velocity.
 * @param t Time in years for the new position
 * @param updateErrors Shall uncertainties be calculated?
 * @param useMaxErrors Use maximum error instead of Gaussian error propagation
 * @return GSL_SUCCESS
 */
int StarGalacticParam::linearOrbit(const mytype t, bool updateErrors, bool useMaxErrors)
{
	x() = x0 + u0*t;
	y() = y0 + v0*t;
	z() = z0 + w0*t;

	if ( updateErrors )
		linearError( t, useMaxErrors );

	time = t;

	return GSL_SUCCESS;
}

/**
 * @brief Calculate errors in position
 * @param t Time in years
 * @param maxError If true, use maximum error instead of Gaussian error propagation
 *
 * Errors in velocities wil not change.
 * @sa StarGalacticParam::epicyclicError()
 */
void StarGalacticParam::linearError(const mytype t, bool maxError)
{
	if ( maxError ) {
		sigmaX = sigmaX0 + t*sigmaU0;
		sigmaY = sigmaY0 + t*sigmaV0;
		sigmaZ = sigmaZ0 + t*sigmaW0;
	} else {
		sigmaX = sqrt( sigmaX0*sigmaX0 + t*sigmaU0*t*sigmaU0 );
		sigmaY = sqrt( sigmaY0*sigmaY0 + t*sigmaV0*t*sigmaV0 );
		sigmaZ = sqrt( sigmaZ0*sigmaZ0 + t*sigmaW0*t*sigmaW0 );
	}
}

/**
 * @brief Calculate position and velocities for a given time
 * @param t Time in years relative to now
 * @return Result of numericOrbit(), epicyclicOrbit(), and linearOrbit(), resp.
 *
 * When using this function always check whether resetODE() should
 * be called first.
 *
 * For the time being this will \b not update the uncertainties!
 */
int StarGalacticParam::orbitAt(const mytype t)
{
	int rc = -1;

	switch ( orbitType ) {
		case Numeric:
			x() = x0; y() = y0; z() = z0;
			u() = u0; v() = v0; w() = w0;
			time = 0;
			rc = numericOrbit( t );
			break;
		case Epicycle:
			rc = epicyclicOrbit( t, false, false );
			break;
		case Linear:
			rc = linearOrbit( t, false, false );
			break;
	}
	return rc;
}

/**
 * @brief Progress the orbit by an amount of time
 * @param t Time to progress in years
 * @return Result of gsl_odeiv2_driver_apply()
 *
 * The function will calculate new coordinates and velocities starting
 * from time point #time to #time+\a t by numerically integrating
 * the equation of motion
 *
 * \sa eqnOfMotion()
 */
int StarGalacticParam::numericOrbit(const mytype t)
{
	int rc = gsl_odeiv2_driver_apply( odeDriver, &time, time+t, cartesian );
//	int rc = gsl_odeiv2_driver_apply_fixed_step( odeDriver, &time, t, 1, cartesian );

	if ( rc != GSL_SUCCESS )
		qWarning( "numericOrbit(): return value = %d", rc );

	return rc;
}

/**
 * @brief Progress the orbit by an amount of time
 * @param t Time to progress in years
 * @return Result of numericOrbit() and epicyclicOrbit(), resp.
 *
 * The function will calculate new coordinates and velocities starting
 * from time point #time to #time+t using a method determined by
 * #orbitType
 *
 * When using this function always check whether resetODE() should
 * be called first.
 */
int StarGalacticParam::progressOrbit(const mytype t)
{
	int rc = -1;

	switch ( orbitType ) {
		case Numeric:
			rc = numericOrbit( t );
			break;
		case Epicycle:
			rc = epicyclicOrbit( time+t, false, false );
			break;
		case Linear:
			rc = linearOrbit( time+t, false, false );
			break;
	}
	return rc;
}

/*!
 * \brief Reset the evolution and stepper objects and the
 * initial integration step size of the GSL ODE driver.
 * \param hstart Step size
 *
 * According to the GSL documentation a reset should be done
 * whenever the next step of the integration will not be the
 * continuation of the previous one.
 *
 * The default value for \a hstart is read from the config file and stored
 * in \a simulParams.StepSize. This function can be used if the integration
 * direction needs to be changed.
 * The value of \a hstart must have the same algebraic sign as the argument
 * of orbitAt() and progressOrbit(), resp.
 *
 * \sa SimulParams
 */
void StarGalacticParam::resetODE(const double hstart)
{
	if ( orbitType != Numeric ) return;

	gsl_odeiv2_driver_reset_hstart( odeDriver, hstart );
}

/**
 * @brief Calculate heliocentric cartesian coordinates
 * @return Matrix with the three components \a x, \a y, \a z
 *         in the first row, and \a u, \a v, \a w in the
 *         second ([pc] and [pc/yr], respectively).
 *
 * For the transformation we first do the translation to
 * the new position of the Sun, and then the coordinate system
 * is rotated around the \a Z axis such that the \a X axis
 * points to the galactic center once again.
 *
 * @sa coordsInertial()
 */
Matrix<mytype> StarGalacticParam::coordsHeliocentric() const
{
	const Vector<mytype>
	        position(3, cartesian),
	        velocity(3, &cartesian[3]);

	Matrix<mytype> result(2,3,0.0);

	switch( orbitType ) {
		case Epicycle:
			// Coordinates are already heliocentric
			result[0] = position;
			result[1] = velocity;
			break;
		case Numeric:
		case Linear:
			// inertial (old) -> rotating (new)

			// Create an instance for the Sun
			StarGalacticParam sun( orbitType, true, time );

			// Now the funny part: we need the base vectors
			// of the new (rotating) system
			// given in inertial system

			// Vector from 'old' origin to galactic center
			// given in inertial system:
			Vector<mytype> oldGalCenter(3);
			oldGalCenter[0] = distGalcenterLSR;

			// Vector from 'old' origin to 'new' origin
			// (its projection into the galactic plane, actually)
			// given in inertial system
			Vector<mytype> newOrigin(3);

			newOrigin[0] = sun.getX();
			newOrigin[1] = sun.getY();
			newOrigin[2] = 0.0;			// Don't change Z direction

			// Vector from 'new' origin to galactic center,
			// given in inertial system:
			Vector<mytype> newGalCenter( oldGalCenter - newOrigin );
#if 0
			const double angle = StarParam::rad2deg( acos(
			        (newGalCenter | oldGalCenter) /
			        (norm2(newGalCenter)*norm2(oldGalCenter)) ) );
			qDebug( "Drehwinkel ist %f Grad", angle );
			qDebug( "Drehwinkel (Epi): %f Grad",
			        StarParam::rad2deg(omega0*time));
#endif
			// Unit vector to galactic center in old (inertial)
			// coordinate system (our new X direction):
			Vector<mytype> e_x( normalize(newGalCenter) );

			// Unit vector of new Z direction (same as the old one)
			// in inertial system:
			Vector<mytype> e_z(3); e_z[2] = 1.0;

			// Unit vector of new Y direction (cross product of e_z and e_x)
			// in inertial system:
			Vector<mytype> e_y( e_z % e_x );

			// Rotation matrix:
			Matrix<mytype> Q(3,3);

			// Filling the rows with the unit vectors will
			// give us a rotation matrix ready for passive
			// rotation (isn't that cute?)
			Q[0] = e_x;
			Q[1] = e_y;
			Q[2] = e_z;

			// For the translation we need
			// the real coordinates of the new origin
			newOrigin[2] = sun.getZ();

			// Now shift and rotate
			result[0] = Q*(position-newOrigin);

			// FIXME: velocity transformation is totally unclear
			Vector<mytype> velSun(3);

			velSun[0] = sun.getU();
			velSun[1] = sun.getV();
			velSun[2] = sun.getW();

			result[1] = Q*(velocity-velSun);

			// TODO: what about uncertainties?
			break;
	}
	return result;
}

/**
 * @brief Get cartesian coordinates relative to the inertial system.
 * @return Matrix with first row containing \a x, \a y, \a z [pc]
 *         and second row \a u, \a v, \a w [pc/yr].
 *
 * @sa coordsHeliocentric()
 */
Matrix<mytype> StarGalacticParam::coordsInertial() const
{
	const Vector<mytype>
	        position(3, cartesian),
	        velocity(3, &cartesian[3]);

	Matrix<mytype> result(2,3,0.0);

	switch( orbitType ) {
		case Numeric:
		case Linear:
			result[0] = position;
			result[1] = velocity;
			break;
		case Epicycle:
#if 1
			// rotating (old) -> inertial (new)

			// Create an instance for the Sun in the
			// inertial system.
			// (Actually, this means we cheat)
			StarGalacticParam sun( Numeric, true, time );

			// Vector from 'new' origin to galactic center
			// in inertial system:
			Vector<mytype> newGalCenter(3);
			newGalCenter[0] = distGalcenterLSR;

			// Vector from 'new' origin to 'old' origin
			// in inertial system:
			Vector<mytype> oldOrigin(3);
			oldOrigin[0] = sun.getX();
			oldOrigin[1] = sun.getY();
			oldOrigin[2] = 0.0;		// Rotation around Z-axis

			// Vector from 'old' origin to galactic center
			// in inertial system:
			Vector<mytype> oldGalCenter( newGalCenter - oldOrigin );

			// Unit vector of 'old' X direction
			// in inertial system:
			Vector<mytype> e_x( normalize( oldGalCenter ) );

			// Unit vector of 'old' Z direction
			// in inertial system:
			Vector<mytype> e_z(3); e_z[2] = 1.0;

			// Unit vector of 'old' Y direction
			// in inertial system
			Vector<mytype> e_y( e_z % e_x );

			// Rotation matrix:
			Matrix<mytype> Q(3,3);

			Q[0] = e_x;		// 1st row
			Q[1] = e_y;		// 2nd row
			Q[2] = e_z;		// 3rd row

			// We need the base vectors as columns
			Q = Q.transpose();

			// The real coordinates of the 'old' origin
			oldOrigin[2] = sun.getZ();

			// Here we go
			result[0] = oldOrigin + (Q*position);

			// The same for the velocity
			Vector<mytype> velSun(3);

			velSun[0] = sun.getU();
			velSun[1] = sun.getV();
			velSun[2] = sun.getW();

			result[1] = velSun + (Q*velocity);
#else
			// The bullshit follows

			// Coordinate system is rotating with angular velocity omega0
			// around galactic center
			const mytype angle = -omega0 * time+0.;
//			qDebug( "Winkel ist %f Grad", StarParam::rad2deg(-angle));

			// Shift origin into galactic center
			location[0] = location[0] - distGalcenterLSR;

			// Rotation matrix
			Matrix<mytype> R(3,3,0.0);
			R[0][0] = cos( angle ); R[0][1] = -sin( angle );
			R[1][0] = -R[0][1];     R[1][1] = R[0][0];
			R[2][2] = 1.0;

			// Rotate
			result[0] = R * location;
			result[1] = R * velocity;

			// Shift back
			result[0][0] += distGalcenterLSR;

			// Add Z component the Sun gained during time
			result[0][2] += kms2pcyr(vLSRz)/nu * sin( nu*(time) ) + lSRz * cos( nu*(time) ) - lSRz;
			result[1][2] += kms2pcyr(vLSRz) * cos( nu*time ) - lSRz * nu * sin( nu*time );

			// Correct for Sun movement
			result[1][0] += kms2pcyr( vLSRx );
			result[1][1] += kms2pcyr( vLSRy + velocityLSROrbitingGalaxy );
			result[1][2] += kms2pcyr( vLSRz );
#endif
			break;
	}

	return result;
}

/**
 * @brief Print the constants used for calculation to standard output
 */
void StarGalacticParam::printConstants(QTextStream &out)
{
	out << ENDL << "Constants used for calculation:" << ENDL << ENDL;

	out << "Coordinate transformation:" << ENDL;
	out << "  distGalcenterLSR          = " << distGalcenterLSR          << ENDL;
	out << "  velocityLSROrbitingGalaxy = " << velocityLSROrbitingGalaxy << ENDL << ENDL;

	out << "LSR:" << ENDL;
	out << QString::fromUtf8( "  X☉ = " ) << lSRx  << ENDL;
	out << QString::fromUtf8( "  Y☉ = " ) << lSRy  << ENDL;
	out << QString::fromUtf8( "  Z☉ = " ) << lSRz  << ENDL;
	out << QString::fromUtf8( "  U☉ = " ) << vLSRx << ENDL;
	out << QString::fromUtf8( "  V☉ = " ) << vLSRy << ENDL;
	out << QString::fromUtf8( "  W☉ = " ) << vLSRz << ENDL << ENDL;

	out << "Epicycles:" << ENDL;
	out << "  omega0 = " << omega0 << ENDL;
	out << "  OortA  = " << OortA  << ENDL;
	out << "  OortB  = " << OortB  << ENDL;
	out << "  kappa  = " << kappa  << ENDL;
	out << "  nu     = " << nu     << ENDL << ENDL;

	out << "Galactic potential:" << ENDL;
	out << "  G       = " <<           galaxyParams.G        << ENDL;
	out << "  Mdisk   = " <<           galaxyParams.Mdisk    << ENDL;
	out << "  Msphere = " <<           galaxyParams.Msphere  << ENDL;
	out << "  a       = " <<           galaxyParams.a        << ENDL;
	out << "  b       = " <<           galaxyParams.b        << ENDL;
	out << "  c       = " <<           galaxyParams.c        << ENDL;
	out << "  d       = " <<           galaxyParams.d        << ENDL;
	out << "  vhalo   = " << pcyr2kms( galaxyParams.vhalo )  << ENDL;
}

void StarGalacticParam::printInfo(QTextStream &out) const
{
	QString oT;
	switch ( getOrbitType() ) {
		case Numeric:
			oT = "Numeric";
			break;
		case Epicycle:
			oT = "Epicycle";
			break;
		case Linear:
			oT = "Linear";
			break;
	}

	out << QString( "%1 (%2):" )
	       .arg( initialValues.getHIP() )
	       .arg( initialValues.getRem() ) << ENDL;
	out << QString( "  Coordinate type: %1" )
	       .arg( initialValues.getCoordType() ) << ENDL;
	out << QString( "  Orbit type     : %1" )
	       .arg( oT ) << ENDL;
}

void StarGalacticParam::printStartValues(QTextStream &out, bool prependHash) const
{
	out << ENDL;
	if ( prependHash ) out << "# ";
	out << "Values from file \""
	    << initialValues.getInputFileName()
	    << "\", line #"
	    << initialValues.getLineNumber() << ":" << ENDL;
	initialValues.printValues( out, prependHash );
}

void StarGalacticParam::printEquatorial(QTextStream &out, bool heliocentric) const
{
	const QString alpha = QString::fromUtf8("α");
	const QString delta = QString::fromUtf8("δ");

	const Vector<mytype> equ( getEquatorial( !heliocentric ) );

	out << QString("  Equatorial coordinates of '%1' at time %L2 yr:")
	       .arg(getHIP())
	       .arg(time, 0, 'f', 0) << ENDL;

	out << QString("    %1 = %2")
	       .arg(alpha)
	       .arg(StarParam::deg2hms( equ[0] ));

	out << QString(";   %1 = %2")
	       .arg(delta)
	       .arg(StarParam::deg2dms( equ[1] ));

	out << QString(";   R = %L1 pc")
	       .arg(equ[2], 0, 'f');

	out << ENDL;
}

/**
 * @brief Calculate galactic coordinates from starting values
 * @return Vector with longitude and latitude in degrees (2 components only!)
 *
 * The function is similar to galacticFromEquatorial() but supposed
 * to be faster because the transformation matrix \b T doesn't
 * have to be initialized each time.
 *
 * See
 * <a href="http://adsabs.harvard.edu/abs/1987AJ.....93..864J">
 * Johnson et al., 1987</a>
 * for details.
 */
Vector<mytype> StarGalacticParam::initLB()
{
	mytype alpha, delta, sinb, cosl, sinl;
	mytype l = INFINITY, b = INFINITY;
	Vector<mytype> equ(3), gal(3), result(2);

	switch ( initialValues.getCoordType() ) {
		case StarParam::Equatorial:
		case StarParam::EquatorialAssoc:
		case StarParam::EquatorialEPM:
		case StarParam::EquatorialEPMAssoc:
			alpha = StarParam::deg2rad( initialValues.getRA() );
			delta = StarParam::deg2rad( initialValues.getDec() );

			equ[0] = cos( delta ) * cos( alpha );
			equ[1] = cos( delta ) * sin( alpha );
			equ[2] = sin( delta );

			gal = T * equ;

			sinb = gal[2];

			// catch rounding errors
			if ( sinb > 1. ) sinb = 1.;
			else if ( sinb < -1. ) sinb = -1.;

			b = asin( sinb );				// b is now in radian

			cosl = gal[0]/* / cos( b )*/;
			sinl = gal[1]/* / cos( b )*/;

			l = StarParam::rad2deg( atan2( sinl, cosl ) );	// in degrees

			if ( l < 0. ) l += 360.;

			b = StarParam::rad2deg( b );	// b is now in degrees
			break;
		case StarParam::GalacticUVW:
		case StarParam::GalacticUVWAssoc:
		case StarParam::GalacticPM:
		case StarParam::GalacticPMAssoc:
			l = initialValues.getLongitude();
			b = initialValues.getLatitude();
			break;
	} /*switch ( initialValues.getCoordType() )*/

	result[0] = l;
	result[1] = b;

	return result;
}

/**
 * @brief Initialize cartesian coordinates
 *
 * Calculate cartesian coordinates from galactic longitude,
 * latitude, and parallax.
 */
void StarGalacticParam::initXYZ()
{
	Vector<mytype> gal( initLB() );

	const mytype l = gal[0];
	const mytype b = gal[1];

	// NOTE: this is not really necessary here
	switch ( initialValues.getCoordType() ) {
		case StarParam::Equatorial:
		case StarParam::GalacticUVW:
		case StarParam::GalacticPM:
		case StarParam::EquatorialEPM:
		case StarParam::EquatorialAssoc:
		case StarParam::GalacticUVWAssoc:
		case StarParam::GalacticPMAssoc:
		case StarParam::EquatorialEPMAssoc:
			const mytype phi = StarParam::deg2rad( l );			// Azimuth angle
			const mytype theta = StarParam::deg2rad( 90. - b );	// Polar angle
			const mytype r = 1000. / initialValues.getParallax();

			x() = r * sin( theta ) * cos( phi );
			y() = r * sin( theta ) * sin( phi );
			z() = r * cos( theta );

			break;
	} /*switch ( initialValues.getCoordType() )*/

#if 0
//	Upper scorpius from Nina
	x() = 140.0;
	y() = -22.0;
	z() =  50.0;
#endif
#if 0
//	Upper scorpius from Bobylev 2008
	x() = 134.0;
	y() = -20.0;
	z() =  52.0;
#endif
}

/**
 * @brief Initialize space velocities
 *
 * Calculate cartesian components of velocities from proper
 * motions and radial velocity. See
 * <a href="http://adsabs.harvard.edu/abs/1987AJ.....93..864J">
 * Johnson et al., 1987</a>
 * for details. Except for epicyclic orbits,
 * the values are then corrected for
 * \li the Solar motion with respect to the Local Standard of Rest
 * @f$v_{lsr}@f$, and
 * \li the Galactic rotational velocity of the Local Standard
 * of Rest @f$v_{gr}@f$.
 *
 * \f$v_{gal} = v_* + v_{lsr} + v_{gr}\f$
 */
void StarGalacticParam::initUVW()
{
	const mytype k = auyr2kms( 1. );
	mytype alpha, delta, theta, phi;
	mytype F_r, F_theta, F_phi;
	Vector<mytype> pm(3), uvw(3);
	Vector<mytype> pmEqu(2), pmEcl(2);

	switch ( initialValues.getCoordType() ) {
		case StarParam::EquatorialEPM:
		case StarParam::EquatorialEPMAssoc:
			alpha = StarParam::deg2rad( initialValues.getRA() );
			delta = StarParam::deg2rad( initialValues.getDec() );

			pmEcl[0] = initialValues.getPMlambda();
			pmEcl[1] = initialValues.getPMbeta();

			pmEqu = equPMFromEclPM(pmEcl, alpha, delta);

			pm[0] = initialValues.getRV();
			pm[1] = k * pmEqu[0] / initialValues.getParallax();
			pm[2] = k * pmEqu[1] / initialValues.getParallax();

			uvw = B( alpha, delta ) * pm;
			break;
		case StarParam::Equatorial:
		case StarParam::EquatorialAssoc:
			alpha = StarParam::deg2rad( initialValues.getRA() );
			delta = StarParam::deg2rad( initialValues.getDec() );

			pm[0] = initialValues.getRV();
			// angles in PMa, PMd, and para must have the same units!
			pm[1] = k * initialValues.getPMa() / initialValues.getParallax();
			pm[2] = k * initialValues.getPMd() / initialValues.getParallax();

			uvw = B( alpha, delta ) * pm;
			break;
		case StarParam::GalacticUVW:
		case StarParam::GalacticUVWAssoc:
			uvw[0] = initialValues.getU();
			uvw[1] = initialValues.getV();
			uvw[2] = initialValues.getW();
			break;
		case StarParam::GalacticPM:
		case StarParam::GalacticPMAssoc:
			theta = StarParam::deg2rad( 90. - initialValues.getLatitude() );
			phi = StarParam::deg2rad( initialValues.getLongitude() );

			F_r = initialValues.getRV();		// in km/s
			F_theta = -k * initialValues.getPMb() / initialValues.getParallax();
			F_phi = k * initialValues.getPMl() / initialValues.getParallax();

			uvw[0] = sin(theta) * cos(phi) * F_r + cos(theta) * cos(phi) * F_theta - sin(phi) * F_phi;
			uvw[1] = sin(theta) * sin(phi) * F_r + cos(theta) * sin(phi) * F_theta + cos(phi) * F_phi;
			uvw[2] = cos(theta) * F_r - sin(theta) * F_theta;
			break;
	} /*switch ( initialValues.getCoordType() )*/

	if ( orbitType != Epicycle ) {
		Vector<mytype>
		        vLSR(3),	// Solar motion with respect to the LSR [km/s]
		        vGR(3);		// Galactic rotational velocity of the LSR [km/s]

		vLSR[0] = vLSRx;	// U component
		vLSR[1] = vLSRy;	// V component
		vLSR[2] = vLSRz;	// W component

		vGR[1] = velocityLSROrbitingGalaxy;

		uvw += vLSR + vGR;
	}

	u() = kms2pcyr( uvw[0] );
	v() = kms2pcyr( uvw[1] );
	w() = kms2pcyr( uvw[2] );
}

/**
 * @brief Calculate uncertainties in \a X, \a Y, \a Z from parallax error
 */
void StarGalacticParam::initSigmaXYZ()
{
	switch ( initialValues.getCoordType() ) {
		case StarParam::Equatorial:
		case StarParam::GalacticUVW:
		case StarParam::GalacticPM:
		case StarParam::EquatorialEPM:
		case StarParam::EquatorialAssoc:
		case StarParam::GalacticUVWAssoc:
		case StarParam::GalacticPMAssoc:
		case StarParam::EquatorialEPMAssoc:
			const mytype pi = initialValues.getParallax();
			const mytype dpi = initialValues.getParaerr();

			sigmaX = fabs( x() * dpi/pi );
			sigmaY = fabs( y() * dpi/pi );
			sigmaZ = fabs( z() * dpi/pi );
			break;
	} /*switch ( initialValues.getCoordType() )*/
}

/**
 * @brief Calculate uncertainties in \a U, \a V, \a W from observation uncertainties
 *
 * See
 * <a href="http://adsabs.harvard.edu/abs/1987AJ.....93..864J">
 * Johnson et al., 1987</a>
 * for details.
 */
void StarGalacticParam::initSigmaUVW()
{
	mytype alpha, delta;
	Matrix<mytype> C;
	Vector<mytype> vec1(3), vec2(3), result(3);
	mytype sigRho, sigPMa, sigPMd, sigPi;
	mytype pi, PMa, PMd;
	const mytype k = auyr2kms( 1. );

	// NOTE: input file gives standard error which we use as sigma

	switch ( initialValues.getCoordType() ) {
		case StarParam::Equatorial:
		case StarParam::EquatorialAssoc:
			alpha = StarParam::deg2rad( initialValues.getRA() );
			delta = StarParam::deg2rad( initialValues.getDec() );

			C = B( alpha,delta );

			sigRho = initialValues.getRVerr();
			sigPMa = initialValues.getPMaerr();
			sigPMd = initialValues.getPMderr();
			sigPi  = initialValues.getParaerr();
			pi     = initialValues.getParallax();
			PMa    = initialValues.getPMa();
			PMd    = initialValues.getPMd();

			vec1[0] = sigRho * sigRho;
			vec1[1] = k*k/pi/pi * ( sigPMa*sigPMa + PMa*PMa*sigPi*sigPi/pi/pi );
			vec1[2] = k*k/pi/pi * ( sigPMd*sigPMd + PMd*PMd*sigPi*sigPi/pi/pi );

			vec2[0] = C[0][1] * C[0][2];
			vec2[1] = C[1][1] * C[1][2];
			vec2[2] = C[2][1] * C[2][2];

			vec2 *= 2*PMa*PMd*k*k*sigPi*sigPi/pi/pi/pi/pi;

			for ( int i=0; i<3; i++ )
				for ( int k=0; k<3; k++ ) {
					C[i][k] = C[i][k] * C[i][k];
				}

			result = C * vec1 + vec2;

			sigmaU = kms2pcyr( sqrt( result[0] ) );
			sigmaV = kms2pcyr( sqrt( result[1] ) );
			sigmaW = kms2pcyr( sqrt( result[2] ) );
			break;
		case StarParam::GalacticUVW:
		case StarParam::GalacticUVWAssoc:
			sigmaU = kms2pcyr( initialValues.getUerr() );
			sigmaV = kms2pcyr( initialValues.getVerr() );
			sigmaW = kms2pcyr( initialValues.getWerr() );
			break;
		case StarParam::GalacticPM:
		case StarParam::GalacticPMAssoc:
		case StarParam::EquatorialEPM:
		case StarParam::EquatorialEPMAssoc:
			// FIXME: initialize sigmas!!!
			sigmaU = sigmaV = sigmaW = 0.;
			break;
	} /*switch ( initialValues.getCoordType() )*/
}

/**
 * @brief
 * Initialize transformation matrix \b T
 * @return 3x3 transformation matrix
 *
 * The matrix is primarily used to transform equatorial
 * coordinates into galactic ones like this:
 *
 * \f$
 * \left[ \begin{array}{c} \!\! \cos b \, \cos l \\ \!\! \cos b \, \sin l \\ \!\! \sin b \end{array} \!\! \right] =
 * \bf{T}
 * \left[ \begin{array}{c} \!\! \cos \delta \, \cos \alpha \\ \!\! \cos \delta \, \sin \alpha \\ \!\! \sin \delta \end{array} \!\! \right]
 * \f$
 *
 * See
 * <a href="http://adsabs.harvard.edu/abs/1987AJ.....93..864J">
 * Johnson et al., 1987</a>
 * for details.
 */
Matrix<mytype> StarGalacticParam::initT()
{
	Matrix<mytype> T1(3,3), T2(3,3), T3(3,3);

	T1[0][0] = cos(theta0);    T1[0][1] =  sin(theta0);   T1[0][2] = 0;
	T1[1][0] = sin(theta0);    T1[1][1] = -cos(theta0);   T1[1][2] = 0;
	T1[2][0] = 0;              T1[2][1] = 0;              T1[2][2] = 1;

	T2[0][0] = -sin(deltaNGP); T2[0][1] =  0;             T2[0][2] = cos(deltaNGP);
	T2[1][0] = 0;              T2[1][1] = -1;             T2[1][2] = 0;
	T2[2][0] = cos(deltaNGP);  T2[2][1] = 0;              T2[2][2] = sin(deltaNGP);

	T3[0][0] = cos(alphaNGP);  T3[0][1] =  sin(alphaNGP); T3[0][2] = 0;
	T3[1][0] = sin(alphaNGP);  T3[1][1] = -cos(alphaNGP); T3[1][2] = 0;
	T3[2][0] = 0;              T3[2][1] = 0;              T3[2][2] = 1;

	return (T1 * T2 * T3);

}

/**
 * @brief Calculate transformation matrix \b B
 * @param alpha Right ascension [rad]
 * @param delta Declination [rad]
 * @return 3x3 transformation matrix
 *
 * The matrix \b B is used to calculate the space-velocity components @a U, @a V, and @a W
 * from proper motion and radial velocity:
 *
 * \f$
 * \left[ \begin{array}{c} \!\! U \\ \!\! V \\ \!\! W \end{array} \!\! \right] =
 * \bf{B} \left[ \begin{array}{c} \!\! \rho \\ \!\! k \mu_\alpha^* / \pi \\ \!\! k \mu_\delta / \pi \end{array} \!\! \right]
 * \f$
 *
 * See
 * <a href="http://adsabs.harvard.edu/abs/1987AJ.....93..864J">
 * Johnson et al., 1987</a>
 * for details.
 */
Matrix<mytype> StarGalacticParam::B(mytype alpha, mytype delta) const
{
	Matrix<mytype> A(3,3);

	const mytype
			sinalpha = sin(alpha),
			cosalpha = cos(alpha),
			sindelta = sin(delta),
			cosdelta = cos(delta);

	A[0][0] = cosalpha*cosdelta;   A[0][1] = -sinalpha;   A[0][2] = -cosalpha * sindelta;
	A[1][0] = sinalpha*cosdelta;   A[1][1] =  cosalpha;   A[1][2] = -sinalpha * sindelta;
	A[2][0] = sindelta;            A[2][1] =  0;          A[2][2] =  cosdelta;

	A = T * A;

	return A;
}

/**
 * @brief Function to calculate the time derivatives of coordinates and velocities
 * according to the equation of motion
 * @param y Array with current coordinates and velocities
 * @param dydt Time derivatives of coordinates and velocites
 * @param params Pointer to parameters defining the gravitational potential
 * @return GSL_SUCCESS or error code
 *
 * The equation of motion
 *
 * \f$
 * \ddot{\vec{r}} = - \vec \nabla \Phi
 * \f$
 *
 * or the system of equations
 *
 * \f$
 * \begin{array}{lll}
 * \dot{x} = u &
 * \dot{y} = v &
 * \dot{z} = w \\
 * \dot{u} = - \displaystyle \frac{\partial \Phi}{\partial x} &
 * \dot{v} = - \displaystyle \frac{\partial \Phi}{\partial y} &
 * \dot{w} = - \displaystyle \frac{\partial \Phi}{\partial z}
 * \end{array}
 * \f$
 *
 * is solved by means of the
 * <a href="http://www.gnu.org/software/gsl/">GNU scientific library</a>.
 *
 * The function is used as a callback and has to be reentrant!
 *
 * The potential that is used here is as follows
 * (see <a href="http://adsabs.harvard.edu/abs/1996ApJ...465..278J"> Johnston et al.</a>
 * for details):
 *
 * \f$
 * \Phi = \Phi_{disk} + \Phi_{sphere} + \Phi_{halo}
 * \f$
 *
 * \f$
 * \begin{array}{ccc}
 * \Phi _{disk} & = & - \displaystyle \frac{GM_{disk}}{\sqrt{x^2 + y^2 + \left(a +\sqrt{z^2 + b^2}\right)^2}} \\ \\
 * \Phi _{sphere} & = & -\displaystyle \frac{GM_{sphere}}{c + \sqrt{x^2 + y^2 + z^2}} \\ \\
 * \Phi _{halo} & = & \displaystyle v_{halo}^2\;\ln (x^2 + y^2 + z^2 +d^2)
 * \end{array}
 * \f$
 *
 * \sa GalaxyParams
 */
int StarGalacticParam::eqnOfMotion(double /*t*/, const double y[], double dydt[], void *params)
{
	const GalaxyParams *p = (GalaxyParams *) params;

	const double
	        X = y[0] - lSRx - distGalcenterLSR,
	        Y = y[1] - lSRy,
	        Z = y[2] - lSRz,
	        U = y[3],
	        V = y[4],
	        W = y[5],

	        A = p->a,
	        B = p->b,
	        C = p->c,
	        D = p->d,
	        G = p->G,
	        Mdisk = p->Mdisk,
	        Msphere = p->Msphere,
	        vhalo = p->vhalo;

	const double R2 = X*X + Y*Y;
	const double r2 = R2 + Z*Z;
	const double r = sqrt( r2 );
	const double beta = sqrt( B*B + Z*Z );

	dydt[0] = U;
	dydt[1] = V;
	dydt[2] = W;

	double denom = (A + beta)*(A + beta) + R2;
	denom *= sqrt( denom );

	const double term1 = G * Mdisk / denom;
	const double term2 = G * Msphere / (r*(C+r)*(C+r));
	const double term3 = vhalo * vhalo * 2.0 / (D*D + r2);

	dydt[3] = -X * (term1 + term2 + term3);
	dydt[4] = -Y * (term1 + term2 + term3);

	denom *= beta;
	const double term4 = G * Mdisk * (A + beta) / denom;

	dydt[5] = -Z * (term4 + term2 + term3);

	return GSL_SUCCESS;
}

/**
 * @brief Jacobian needed by the
 *        <a href="http://www.gnu.org/software/gsl/">GSL</a>
 *        for solving the differential equation
 * @param y Current positions and velocities
 * @param dfdy Where to store the jacobian
 * @param dfdt Where to store the time derivatives
 * @param params Pointer to parameters defining the gravitational potential
 * @return GSL_SUCCESS or error code
 *
 * Note that Runge-Kutta-Fehlberg (4,5) does \b not use a Jacobian!
 *
 * The function is used as a callback and has to be reentrant!
 * \sa eqnOfMotion()
 *
 * Doxygen bug? Ignore the following reference to b!
 */
int StarGalacticParam::jacEqnOfMotion(double /*t*/, const double y[], double *dfdy, double dfdt[], void *params)
{
	const GalaxyParams *p = (GalaxyParams *) params;

	const double
	        X = y[0] - distGalcenterLSR,
	        Y = y[1],
	        Z = y[2],

	        a = p->a,
	        b = p->b,
	        c = p->c,
	        d = p->d,
	        G = p->G,
	        Mdisk = p->Mdisk,
	        Msphere = p->Msphere,
	        vhalo = p->vhalo;

	const double R2 = X*X + Y*Y;
	const double r2 = R2 + Z*Z;
	const double r = sqrt( r2 );
	const double beta2 = b*b + Z*Z;
	const double beta = sqrt( beta2 );

	gsl_matrix_view dfdy_mat
	        = gsl_matrix_view_array (dfdy, 6, 6);
	gsl_matrix *m = &dfdy_mat.matrix;

	gsl_matrix_set_zero( m );

	gsl_matrix_set( m, 0, 3, 1.0 );
	gsl_matrix_set( m, 1, 4, 1.0 );
	gsl_matrix_set( m, 2, 5, 1.0 );

//	J_41
	double denom1 = (a+beta)*(a+beta) + R2;
	double num1 = G * Mdisk* (denom1 - X*X*3.0);
	denom1 *= denom1 * sqrt( denom1 );
	double term1 = num1 / denom1;

	double denom2 = r2 * r * (c+r)*(c+r)*(c+r);
	double term2 = (G * Msphere * (r2 * (c+r) - (c+r*3.0) * X*X)) / denom2;

	double denom3 = d*d + r2;
	double num3 = vhalo*vhalo * 2.0 * (denom3 - X*X*2.0);
	denom3 *= denom3;

	double term3 = num3 / denom3;

	gsl_matrix_set( m, 3, 0, -term1 - term2 - term3 );

//	J_42, J_51
	term1 = G * Mdisk * 3.0 / denom1;
	term2 = G * Msphere * (c + r*3.0) / denom2;
	term3 = vhalo*vhalo*4.0 / denom3;

	double tmp = X * Y * (term1 + term2 + term3);
	gsl_matrix_set( m, 3, 1, tmp );
	gsl_matrix_set( m, 4, 0, tmp );

//	J_43, J_61, J_53, J_62
	term1 *= ((a + beta)/beta);

	tmp = Z * (term1 + term2 + term3);
	gsl_matrix_set( m, 3, 2, X * tmp );
	gsl_matrix_set( m, 5, 0, X * tmp );
	gsl_matrix_set( m, 4, 2, Y * tmp );
	gsl_matrix_set( m, 5, 1, Y * tmp );

//	J_52
	denom1 = (a+beta)*(a+beta) + R2;
	num1 = G * Mdisk * (denom1 - Y*Y*3.0);
	denom1 *= denom1 * sqrt( denom1 );
	term1 = num1 / denom1;
	term2 = (G * Msphere * (r2 * (c+r) - (c+r*3.0) * Y*Y)) / denom2;
	term3 = vhalo*vhalo * 2.0 * (d*d + r2 - Y*Y*2.0) / denom3;

	gsl_matrix_set( m, 4, 1, -term1 - term2 - term3 );

//	J_63
	const double beta3 = beta2 * beta;
	const double beta4 = beta2 * beta2;
	const double beta5 = beta4 * beta;
	const double a2 = a*a;
	const double Z2 = Z*Z;

	denom1 *= beta3;
	num1 = G*Mdisk*(-a*(a2+R2)*Z2 - a2*Z2*beta*5.0 + a*(a2+R2-Z2*7.0)*beta2 + (a2*3.0+R2-Z2*3.0)*beta3 + a*beta4*3.0 + beta5);
	term1 = num1 / denom1;

	const double num2 = G * Msphere * (r2*(c+r) - (c+r*3.0)*Z2);
	term2 = num2 / denom2;

	num3 = vhalo * vhalo * 2.0 * (d*d + r2 - Z2*2.0);
	term3 = num3 / denom3;

	gsl_matrix_set( m, 5, 2, -term1 - term2 - term3 );

//	dfdt
	for ( int i = 0; i < 6; i++ )
		dfdt[i] = 0.0;

	return GSL_SUCCESS;
}

/*!
 * \brief Initialize #odeSys and #odeDriver
 */
void StarGalacticParam::initODE()
{
	odeSys.function = eqnOfMotion;
	odeSys.jacobian = jacEqnOfMotion;	// But Runge-Kutta-45 doesn't use it
	odeSys.dimension = 6;
	odeSys.params = &galaxyParams;
	odeDriver = gsl_odeiv2_driver_alloc_y_new( &odeSys, ODESTEPPER, simulParams.StepSize, 1e-6, 0.0 );
}

void StarGalacticParam::printGalactic(QTextStream &out, bool heliocentric) const
{
	const QString deg = QString::fromUtf8("°");

	const Vector<mytype> gal( getGalactic( !heliocentric ) );

	out << QString("  Galactic coordinates of '%1' at time %L2 yr:")
	       .arg(getHIP())
	       .arg(time, 0, 'f', 0) << ENDL;

	out << QString("    l = %L1%2")
	       .arg((double) gal[0], 0, 'f')
	       .arg(deg);

	out << QString(";   b = %L1%2")
	       .arg((double) gal[1], 0, 'f')
	       .arg(deg);

	out << QString(";   R = %L1 pc")
	       .arg((double) gal[2], 0, 'f');

	out << ENDL;
}

void StarGalacticParam::printXYZ(QTextStream &out, bool heliocentric) const
{
	const QString pm = QString::fromUtf8( "±" );
	const Vector<mytype> xyz( heliocentric ? coordsHeliocentric()[0] : coordsInertial()[0] );

	out << QString("  %3Cartesian coordinates of '%1' at time %L2 yr:")
	       .arg(getHIP())
	       .arg(time, 0, 'f', 0)
	       .arg(heliocentric ? "Heliocentric " : "") << ENDL;

	out << QString("    X = (%L1 %2 %L3) pc")
	       .arg((double) xyz[0])
	       .arg(pm)
	       .arg((double) sigmaX) << ENDL;
	out << QString("    Y = (%L1 %2 %L3) pc")
	       .arg((double) xyz[1])
	       .arg(pm)
	       .arg((double) sigmaY) << ENDL;
	out << QString("    Z = (%L1 %2 %L3) pc")
	       .arg((double) xyz[2])
	       .arg(pm)
	       .arg((double) sigmaZ) << ENDL;
}

void StarGalacticParam::printUVW(QTextStream &out, bool heliocentric) const
{
	const QString pm = QString::fromUtf8( "±" );
	const Vector<mytype> velocity( heliocentric ? coordsHeliocentric()[1] : coordsInertial()[1] );

	const mytype
	        u = pcyr2kms(velocity[0]),
	        v = pcyr2kms(velocity[1]),
	        w = pcyr2kms(velocity[2]);

	out << QString("  Space velocity components %1of '%2' at time %L3 yr:")
	       .arg( (heliocentric||orbitType==Epicycle) ? "(relative to Sun) " : "" )
	       .arg(getHIP())
	       .arg(time, 0, 'f', 0) << ENDL;

	out << QString("    U = (%L1 %2 %L3) km/s")
	       .arg((double) u )
	       .arg(pm)
	       .arg((double) pcyr2kms( sigmaU )) << ENDL;

	out << QString("    V = (%L1 %2 %L3) km/s")
	       .arg((double) v )
	       .arg(pm)
	       .arg((double) pcyr2kms( sigmaV )) << ENDL;

	out << QString("    W = (%L1 %2 %L3) km/s")
	       .arg((double) w )
	       .arg(pm)
	       .arg((double) pcyr2kms( sigmaW )) << ENDL;
}

/**
 * @brief Formatted output of proper motion and radial velocity
 * @param out Text stream for output
 * @param heliocentric Relative to Sun or coordinate origin.
 *        Setting this to \b false is adventurous.
 */
void StarGalacticParam::printPM(QTextStream &out, bool heliocentric) const
{
	const QString mu = QString::fromUtf8( "μ" );
	const QString alpha = QString::fromUtf8( "α" );
	const QString delta = QString::fromUtf8( "δ" );
	const QString rho = QString::fromUtf8( "ρ" );

	const Matrix<mytype> coords( heliocentric ? coordsHeliocentric() : coordsInertial() );

	const Vector<mytype> gal( lbRFromXYZ( coords[0] ) );
	const Vector<mytype> equ( equatorialFromGalactic( gal ) );

	const Vector<mytype> pm(
	            pmFromVelocities( StarParam::deg2rad( equ[0] ),
	                              StarParam::deg2rad( equ[1] ),
	                              1./equ[2],
	                              pcyr2kms( coords[1][0] ),
	                              pcyr2kms( coords[1][1] ),
	                              pcyr2kms( coords[1][2] ) )
	                       );

	out << QString("  Proper motion of '%1' at time %L2 yr:")
	       .arg(getHIP())
	       .arg(time, 0, 'f', 0) << ENDL;
	out << QString("    %1    = %L2 km/s").arg(rho).arg(pm[0]) << ENDL;
	out << QString("    %1_%2* = %L3 mas/yr").arg(mu).arg(alpha).arg(pm[1]*1000.) << ENDL;
	out << QString("    %1_%2  = %L3 mas/yr").arg(mu).arg(delta).arg(pm[2]*1000.) << ENDL;
}

void StarGalacticParam::printSigmas(QTextStream &out) const
{
	const QString sg = QString::fromUtf8( "σ" );

	out << QString("Location errors of '%1' at time %L2 yr:")
	       .arg(getHIP())
	       .arg(time, 0, 'f', 0) << ENDL;
	out << QString("  %1x0 = %L2 pc").arg(sg).arg(sigmaX0);
	out << QString("  %1y0 = %L2 pc").arg(sg).arg(sigmaY0);
	out << QString("  %1z0 = %L2 pc").arg(sg).arg(sigmaZ0) << ENDL;
	out << QString("  %1x  = %L2 pc").arg(sg).arg(sigmaX);
	out << QString("  %1y  = %L2 pc").arg(sg).arg(sigmaY);
	out << QString("  %1z  = %L2 pc").arg(sg).arg(sigmaZ) << ENDL;
	out << "Velocity errors:" << ENDL;
	out << QString("  %1u0 = %L2 km/s").arg(sg).arg(pcyr2kms(sigmaU0));
	out << QString("  %1v0 = %L2 km/s").arg(sg).arg(pcyr2kms(sigmaV0));
	out << QString("  %1w0 = %L2 km/s").arg(sg).arg(pcyr2kms(sigmaW0)) << ENDL;
	out << QString("  %1u  = %L2 km/s").arg(sg).arg(pcyr2kms(sigmaU));
	out << QString("  %1v  = %L2 km/s").arg(sg).arg(pcyr2kms(sigmaV));
	out << QString("  %1w  = %L2 km/s").arg(sg).arg(pcyr2kms(sigmaW)) << ENDL;
}

/**
 * @brief Three dimensional distance between two stars
 * @param st Parameters of second star
 * @return Distance in parsec
 *
 * Calculates distance according to
 *
 * \f$ \displaystyle
 * d = \sqrt{(x_2-x_1)^2+(y_2-y_1)^2+(z_2-z_1)^2}
 * \f$
 *
 * To be used like:
 * @code StarGalacticParam star1, star2;
 *
 * mytype dist = star1||star2 @endcode
 */
mytype StarGalacticParam::operator ||(const StarGalacticParam &st) const
{
	mytype dx = getX() - st.getX();
	mytype dy = getY() - st.getY();
	mytype dz = getZ() - st.getZ();

	return sqrt( dx*dx + dy*dy + dz*dz );
}

/**
 * @brief Distance to another star in one dimension
 * @param st Parameters of second star
 * @param direction Either 'X', 'Y', or 'Z'
 * @return Distance in parsec
 */
mytype StarGalacticParam::distanceTo(const StarGalacticParam &st, char direction) const
{
	mytype dist;

	switch( direction ) {
		case 'x':
		case 'X': dist = fabs( st.getX() - getX() ); break;
		case 'y':
		case 'Y': dist = fabs( st.getY() - getY() ); break;
		case 'z':
		case 'Z': dist = fabs( st.getZ() - getZ() ); break;
		default:  dist = 0;                break;
	}

	return dist;
}

/**
 * @brief Check if we are close to another star
 * @param st Second star
 * @param withinSigma Confidence interval (e.g. 3 for 3σ)
 * @return @a True or @a false
 */
bool StarGalacticParam::isCloseTo(const StarGalacticParam &st, double withinSigma) const
{
	return  distanceTo( st, 'X') <= withinSigma*(sigmaX + st.sigmaX) &&
			distanceTo( st, 'Y') <= withinSigma*(sigmaY + st.sigmaY) &&
			distanceTo( st, 'Z') <= withinSigma*(sigmaZ + st.sigmaZ);
}


void StarGalacticParam::writeConfig(QString fName)
{
	QSettings settings( fName, QSettings::IniFormat );
	settings.beginGroup( "Coordinate_Transformation");
	settings.setValue( "distGalcenterLSR", distGalcenterLSR );
	settings.setValue( "velocityLSROrbitingGalaxy", velocityLSROrbitingGalaxy );
	settings.endGroup();

	settings.beginGroup( "LSR" );
	settings.setValue( "x", lSRx );
	settings.setValue( "y", lSRy );
	settings.setValue( "z", lSRz );
	settings.setValue( "v_x", vLSRx );
	settings.setValue( "v_y", vLSRy );
	settings.setValue( "v_z", vLSRz );
	settings.endGroup();

	settings.beginGroup( "Epicycle_Constants" );
	settings.setValue( "omega0", omega0 );
	settings.setValue( "OortA", OortA );
	settings.setValue( "OortB", OortB );
	settings.setValue( "kappa", kappa );
	settings.setValue( "nu", nu );
	settings.endGroup();

	settings.beginGroup( "Galactic_Potential");
	settings.setValue( "G", galaxyParams.G );
	settings.setValue( "Mdisk", galaxyParams.Mdisk );
	settings.setValue( "Msphere", galaxyParams.Msphere );
	settings.setValue( "a", galaxyParams.a );
	settings.setValue( "b", galaxyParams.b );
	settings.setValue( "c", galaxyParams.c );
	settings.setValue( "d", galaxyParams.d );
	settings.setValue( "vhalo", StarGalacticParam::pcyr2kms( galaxyParams.vhalo ) );
	settings.endGroup();

	settings.beginGroup( "Simulation" );
	settings.setValue( "Threads",       simulParams.Threads );
	settings.setValue( "Orbits",        simulParams.Orbits );
	settings.setValue( "Steps",         simulParams.Steps );
	settings.setValue( "StepSize",      simulParams.StepSize);
	settings.setValue( "Width",         simulParams.Width );
	settings.setValue( "Star1",         simulParams.Star1 );
	settings.setValue( "Star2",         simulParams.Star2 );
	settings.setValue( "Assoc",         simulParams.Assoc );
	settings.setValue( "Type",          simulParams.Type );
	settings.setValue( "Criterion",     simulParams.Criterion );
	settings.setValue( "Limit",         simulParams.Limit );
	settings.setValue( "InFile",        simulParams.InFile );
	settings.setValue( "OutFile",       simulParams.OutFile );
	settings.setValue( "RemoveIfEmpty", simulParams.RemoveFile ? "Yes" : "No" );
	settings.endGroup();
}


void StarGalacticParam::readConfig(QString fName)
{
	if ( !QFile::exists( fName ) )
		qFatal( "Error: config file '%s' not found!\n", fName.toLatin1().data() );

	QSettings settings( fName, QSettings::IniFormat );

	settings.beginGroup( "Coordinate_Transformation");
	distGalcenterLSR          = settings.value( "distGalcenterLSR", 0 ).toDouble();
	velocityLSROrbitingGalaxy = settings.value( "velocityLSROrbitingGalaxy", 0 ).toDouble();
	settings.endGroup();

	settings.beginGroup( "LSR" );
	lSRx  = settings.value( "x", 0 ).toDouble();
	lSRy  = settings.value( "y", 0 ).toDouble();
	lSRz  = settings.value( "z", 0 ).toDouble();
	vLSRx = settings.value( "v_x", 0 ).toDouble();
	vLSRy = settings.value( "v_y", 0 ).toDouble();
	vLSRz = settings.value( "v_z", 0 ).toDouble();
	settings.endGroup();

	settings.beginGroup( "Epicycle_Constants" );
	omega0 = settings.value( "omega0", 0 ).toDouble();
	OortA  = settings.value( "OortA", 0 ).toDouble();
	OortB  = settings.value( "OortB", 0 ).toDouble();
	kappa  = settings.value( "kappa", 0 ).toDouble();
	nu     = settings.value( "nu", 0 ).toDouble();
	settings.endGroup();

	settings.beginGroup( "Galactic_Potential");
	galaxyParams.G       = settings.value( "G", 0 ).toDouble();
	galaxyParams.Mdisk   = settings.value( "Mdisk", 0 ).toDouble();
	galaxyParams.Msphere = settings.value( "Msphere", 0 ).toDouble();
	galaxyParams.a       = settings.value( "a", 0 ).toDouble();
	galaxyParams.b       = settings.value( "b", 0 ).toDouble();
	galaxyParams.c       = settings.value( "c", 0 ).toDouble();
	galaxyParams.d       = settings.value( "d", 0 ).toDouble();
	galaxyParams.vhalo   = StarGalacticParam::kms2pcyr( settings.value( "vhalo", 0 ).toDouble() );
	settings.endGroup();

	settings.beginGroup( "Simulation" );
	simulParams.Threads    = settings.value( "Threads", 0 ).toInt();
	simulParams.Orbits     = settings.value( "Orbits", 1000 ).toInt();
	simulParams.Steps      = settings.value( "Steps", 1000 ).toInt();
	simulParams.Width      = settings.value( "Width", 10000 ).toInt();
	simulParams.StepSize   = settings.value( "StepSize" ).toDouble();
	simulParams.Type       = settings.value( "Type", "Potential" ).toString();
	simulParams.Criterion  = settings.value( "Criterion" ).toString();
	simulParams.Limit      = settings.value( "Limit", 10.0 ).toDouble();
	simulParams.Star1      = settings.value( "Star1" ).toString();
	simulParams.Star2      = settings.value( "Star2" ).toString();
	simulParams.Assoc      = settings.value( "Assoc" ).toString();
	simulParams.InFile     = settings.value( "InFile", "input.txt" ).toString();
	simulParams.OutFile    = settings.value( "OutFile", "auto" ).toString();
	simulParams.RemoveFile = settings.value( "RemoveIfEmpty", "No").toString().toUpper() == "YES";
	settings.endGroup();

	if ( (simulParams.Criterion != CRITTETZLAFF) &&
	     (simulParams.Criterion != CRITHOOGERWERF) &&
	     (simulParams.Criterion != CRITRALPH1) &&
	     (simulParams.Criterion != CRITVARLERI1) ) {
		qDebug( "Error in config file: Criterion must be either \"Tetzlaff\" or \"Hoogerwerf\" or \"Ralph1\" or \"Valeri1\"\n");
		exit( 2 );
	}

	if ( (simulParams.Assoc == 0) && (simulParams.Criterion == CRITTETZLAFF) ) {
		qDebug( "Error in config file: \"Tetzlaff\" criterion needs association\n" );
		exit( 2 );
	}
}
