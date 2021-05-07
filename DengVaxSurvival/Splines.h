#ifndef SPLINES_HEADER_INCLUDED
#define SPLINES_HEADER_INCLUDED

#include "StructureDefs.h"

DType Hill						(DType TimePostDose, DType HillCoeff, DType Halflife);
void CoefficientsFromKnots		(DType* xKnots, DType* yKnots, int FirstKnot, DType *Coeffs);
void CoefficientsFromKnots_Line	(DType* xKnots, DType* yKnots, int FirstKnot, DType *Coeffs);
void CoefficientsFromKnots_Cubic(DType* xKnots, DType* yKnots, int FirstKnot, DType *Coeffs);
int FindLastPolyBefore_t		(DType t, DType *xKnots, const int &PolynomialsPerCountry);
inline DType Spline				(DType xValue, DType* xKnots, DType* yKnots, DType **Coeffs, int KnotsPerSpline, int PolynomialsPerSpline);
inline DType Line				(DType xValue, DType* xKnots, DType* yKnots, DType **Coeffs, int KnotsPerSpline, int PolynomialsPerSpline);
inline DType Step				(DType xValue, DType* xKnots, DType* yKnots, DType **Coeffs, int KnotsPerSpline, int PolynomialsPerSpline);
inline DType Cubic				(DType xValue, DType* xKnots, DType* yKnots, DType **Coeffs, int KnotsPerSpline, int PolynomialsPerSpline);
DType BaselineHazard			(DType t, const int &country, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType WaningDurationSpline		(int AgeInYears, int BaselineSeroStatus, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType AS_HazMultSpline			(int AgeInYears, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);

DType EfficacyMultiplierSpline	(DType age, int BaselineSeroStatus, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType IntBaseHaz				(int Lower_t, int Upper_t, const int &country, const Params_Struct &PARAMS);
DType IntVacHaz					(int Lower_t, int Upper_t, int Dose_1_Day, int Dose_2_Day, int Dose_3_Da, int AgeInYears, char BaselineSeroStatus, const int &country, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType VaccineHazard				(int CalendarDay, int Dose_1_Day, int Dose_2_Day, int Dose_3_Day, int AgeInYears, char BaselineSeroStatus, const int &country, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);
DType IntVacHaz_SingleDay		(int CalendarDay, int Dose_1_Day, int Dose_2_Day, int Dose_3_Day, int AgeInYears, char BaselineSeroStatus, const int &country, const Params_Struct &PARAMS, const Housekeeping_Struct &HOUSE);

#define _I64(x) static_cast<int64_t>(x)

#endif