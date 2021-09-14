#ifndef ACS_DATAGEN_H
#define ACS_DATAGEN_H
#include "macros.h"
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include "bessel.h"

/**
 * @brief This variable is unset by the ACS thread at first execution.
 * 
 */
volatile int first_run = 1;

/**
 * @brief Creates buffer for \f$\vec{\omega}\f$.
 * 
 */
DECLARE_BUFFER(g_W, float); // omega global circular buffer
/**
 * @brief Creates buffer for \f$\vec{B}\f$.
 * 
 */
DECLARE_BUFFER(g_B, double); // magnetic field global circular buffer
/**
 * @brief Creates buffer for \f$\vec{\dot{B}}\f$.
 * 
 */
DECLARE_BUFFER(g_Bt, double); // Bdot global circular buffer1
/**
 * @brief Creates vector for target angular momentum.
 * 
 */
DECLARE_VECTOR(g_L_target, float); // angular momentum target vector
/**
 * @brief Creates vector for target angular speed.
 * 
 */
DECLARE_VECTOR(g_W_target, float); // angular velocity target vector
/**
 * @brief Creates buffer for sun vector.
 * 
 */
DECLARE_BUFFER(g_S, float); // sun vector
/**
 * @brief Storage for current coarse sun sensor lux measurements.
 * 
 */
float g_CSS[7]; // current CSS lux values, in HITL this will be populated by TSL2561 code
/**
 * @brief Indicate if Mux channel has error
 * 
 */
bool mux_err_channel[3] = {false, false, false};
/**
 * @brief Storage for current fine sun sensor angle measurements.
 * 
 */
float g_FSS[2]; // current FSS angles, in rad; in HITL this will be populated by NANOSSOC A60 driver
/**
 * @brief Stores the return value of FSS algorithm
 * 
 */
int g_FSS_RET; // return value of FSS
/**
 * @brief Current index of the \f$\vec{B}\f$ circular buffer.
 * 
 */
int mag_index = -1;
/**
 * @brief Current index of the \f$\vec{\omega}\f$ circular buffer.
 * 
 */
int omega_index = -1;
/**
 * @brief Current index of the \f$\vec{\dot{B}}\f$ circular buffer.
 * 
 */
int bdot_index = -1;
/**
 * @brief Current index of the sun vector circular buffer.
 * 
 */
int sol_index = -1; // circular buffer indices, -1 indicates uninitiated buffer
/**
 * @brief Indicates if the \f$\vec{B}\f$ circular buffer is full.
 * 
 */
int B_full = 0;
/**
 * @brief Indicates if the \f$\vec{\dot{B}}\f$ circular buffer is full.
 * 
 */
int Bdot_full = 0;
/**
 * @brief Indicates if the \f$\vec{\omega}\f$ circular buffer is full.
 * 
 */
int W_full = 0;
/**
 * @brief Indicates if the sun vector circular buffer is full.
 * 
 */
int S_full = 0; // required to deal with the circular buffer problem
/**
 * @brief This variable is set by checkTransition() if the satellite does not detect the sun.
 * 
 */
uint8_t g_night = 0; // night mode?
/**
 * @brief This variable contains the current state of the flight system.
 * 
 */
uint8_t g_acs_mode = 0; // Detumble by default
/**
 * @brief This variable is unset when the system is detumbled for the first time after a power cycle.
 * 
 */
uint8_t g_first_detumble = 1; // first time detumble by default even at night
/**
 * @brief Counts the number of cycles on the ACS thread.
 * 
 */
unsigned long long acs_ct = 0; // counts the number of ACS steps
/**
 * @brief Moment of inertia of the satellite (SI).
 * 
 */
float MOI[3][3] = {{0.0821, 0, 0},
                   {0, 0.0752, 0},
                   {0, 0, 0.0874}};
/**
 * @brief Inverse of the moment of inertia of the satellite (SI).
 * 
 */
float IMOI[3][3] = {{12.1733, 0, 0},
                    {0, 13.2941, 0},
                    {0, 0, 11.4661}};
/**
 * @brief Current timestamp after readSensors() in ACS thread, used to keep track of time taken by ACS loop.
 * 
 */
unsigned long long g_t_acs;

/**
 * @brief Dipole moment of the magnetorquer rods
 * 
 */
static float DIPOLE_MOMENT = 0.22; // A m^-2
/**
 * @brief ACS loop time period
 * 
 */
static uint32_t DETUMBLE_TIME_STEP = 100000; // 100 ms for full loop

void getOmega(void)
{
    if (mag_index < 2 && B_full == 0) // not enough measurements
        return;
    // once we have measurements, we declare that we proceed
    if (omega_index == SH_BUFFER_SIZE - 1) // hit max, buffer full
        W_full = 1;
    omega_index = (1 + omega_index) % SH_BUFFER_SIZE;                             // calculate new index in the circular buffer
    int8_t m0, m1;                                                                // temporary addresses
    m1 = bdot_index;                                                              // current address
    m0 = (bdot_index - 1) < 0 ? SH_BUFFER_SIZE - bdot_index - 1 : bdot_index - 1; // previous address, wrapped around the circular buffer
    float freq;
    freq = 1e6 / DETUMBLE_TIME_STEP;                     // time units!
    CROSS_PRODUCT(g_W[omega_index], g_Bt[m1], g_Bt[m0]); // apply cross product
    float norm2 = NORM2(g_Bt[m0]);
    VECTOR_MIXED(g_W[omega_index], g_W[omega_index], freq / norm2, *); // omega = (B_t dot x B_t-dt dot)*freq/Norm2(B_t dot)
    // Apply correction // There is fast runaway with this on
    // DECLARE_VECTOR(omega_corr0, float);                            // declare temporary space for correction vector
    // MATVECMUL(omega_corr0, MOI, g_W[m1]);                          // MOI X w[t-1]
    // DECLARE_VECTOR(omega_corr1, float);                            // declare temporary space for correction vector
    // CROSS_PRODUCT(omega_corr1, g_W[m1], omega_corr0);              // store into temp 1
    // MATVECMUL(omega_corr1, IMOI, omega_corr0);                     // store back into temp 0
    // VECTOR_MIXED(omega_corr1, omega_corr1, -freq, *);              // omega_corr = freq*(MOI-1)*(-w[t-1] X MOI*w[t-1])
    // VECTOR_OP(g_W[omega_index], g_W[omega_index], omega_corr1, +); // add the correction term to omega
    APPLY_FBESSEL(g_W, omega_index); // Bessel filter of order 3
    return;
}

void getSVec(void)
{
    if (sol_index == SH_BUFFER_SIZE - 1) // hit max, buffer full
        S_full = 1;
    sol_index = (sol_index + 1) % SH_BUFFER_SIZE;
#ifndef M_PI
/**
 * @brief Approximate definition of Pi in case M_PI is not included from math.h
 */
#define M_PI 3.1415
#endif
    // check if FSS results are acceptable
    // if they are, use that to calculate the sun vector
    // printf("[FSS] %.3f %.3f\n", fsx * 180. / M_PI, fsy * 180. / M_PI);

    // get average -Z luminosity from 2 sensors
    float znavg = 0;
    for (int i = 5; i < 7; i++)
        znavg += g_CSS[i];
    znavg *= 0.5f;

    x_g_S[sol_index] = g_CSS[0] - g_CSS[1]; // +x - -x
    y_g_S[sol_index] = g_CSS[2] - g_CSS[3]; // +x - -x
    z_g_S[sol_index] = g_CSS[4] - znavg;    // +z - avg(-z)

    float css_mag = NORM(g_S[sol_index]); // norm of the CSS lux values
#define CSS_MIN_LUX_THRESHOLD 500
    if (css_mag < CSS_MIN_LUX_THRESHOLD) // night time logic
    {
        g_night = 1;
        VECTOR_CLEAR(g_S[sol_index]); // return 0 solar vector
#ifdef ACS_PRINT
        printf("[" RED "FSS" RST "]");
#endif // ACS_PRINT
    }
    else
    {
        g_night = 0;
        NORMALIZE(g_S[sol_index], g_S[sol_index]); // return normalized sun vector
#ifdef ACS_PRINT
        printf("[" YLW "FSS" RST "]");
#endif // ACS_PRINT
    }
    printf("[sunvec %d] %0.3f %0.3f %0.3f\n", sol_index, x_g_S[sol_index], y_g_S[sol_index], z_g_S[sol_index]);
    return;
}

int readSensors(void)
{
    // read magfield, CSS, FSS
    static double tnow = 0;
    if (mag_index >= 0)
        printf("In readSensors(): acs count %llu, mag_index %d, Bx %lf By %lf Bz %lf tnow %lf...\n", acs_ct, mag_index, x_g_B[mag_index], y_g_B[mag_index], z_g_B[mag_index], tnow);
    acs_ct++;
    tnow += 0.1; // 0.1 seconds
    int status = 1;
    if (mag_index == SH_BUFFER_SIZE - 1) // hit max, buffer full
        B_full = 1;
    mag_index = (mag_index + 1) % SH_BUFFER_SIZE;
    VECTOR_CLEAR(g_B[mag_index]); // clear the current B                                                  /
    // HITL
    double mag_measure[3];
    mag_measure[0] = 50 * sin(tnow * 0.5) + (((1.0 * rand()) / RAND_MAX) - 0.5); // noise
    mag_measure[1] = 50 * cos(tnow * 0.5) + (((1.0 * rand()) / RAND_MAX) - 0.5); // noise
    mag_measure[2] = ((1.0 * rand()) / RAND_MAX - 0.5); // noise + sine
    // read coarse sun sensors
    double sun_ang = sin(tnow * 0.1) * 15 + 30;
    g_CSS[0] = 7000 * sin(sun_ang * M_PI / 180) * cos(tnow * 0.5) + ((100.0 * rand()) / RAND_MAX - 50);
    g_CSS[1] = -g_CSS[0];
    g_CSS[2] = 7000 * sin(sun_ang * M_PI / 180) * sin(tnow * 0.5) + ((100.0 * rand()) / RAND_MAX - 50);
    g_CSS[3] = -g_CSS[2];
    g_CSS[4] = 7000 * cos(sun_ang * M_PI / 180) + ((100.0 * rand()) / RAND_MAX - 50);
    g_CSS[5] = -g_CSS[4];
    g_CSS[6] = g_CSS[5];

    DECLARE_VECTOR(mag_mes, double);
    x_mag_mes = mag_measure[0]; // / 6.842;
    y_mag_mes = mag_measure[1]; // / 6.842;
    z_mag_mes = mag_measure[2]; // / 6.842;
    DECLARE_VECTOR(mag_val, double);
    NORMALIZE(mag_val, mag_mes);
    VECTOR_MIXED(mag_val, mag_val, 600, *);
    x_g_B[mag_index] = x_mag_mes; // scaled to milliGauss
    y_g_B[mag_index] = y_mag_mes;
    z_g_B[mag_index] = z_mag_mes;
    APPLY_DBESSEL(g_B, mag_index); // bessel filter

    // printf("readSensors: Bx: %f By: %f Bz: %f\n", x_g_B[mag_index], y_g_B[mag_index], z_g_B[mag_index]);
    // put values into g_Bx, g_By and g_Bz at [mag_index] and takes 18 ms to do so (implemented using sleep)
    if (mag_index < 1 && B_full == 0)
        return status;
    // if we have > 1 values, calculate Bdot
    if (bdot_index == SH_BUFFER_SIZE - 1) // hit max, buffer full
        B_full = 1;
    bdot_index = (bdot_index + 1) % SH_BUFFER_SIZE;
    int8_t m0, m1;
    m1 = mag_index;
    m0 = (mag_index - 1) < 0 ? SH_BUFFER_SIZE - mag_index - 1 : mag_index - 1;
    double freq = 1e6 / (DETUMBLE_TIME_STEP * 1.0);
    VECTOR_OP(g_Bt[bdot_index], g_B[m1], g_B[m0], -);
    VECTOR_MIXED(g_Bt[bdot_index], g_Bt[bdot_index], freq, *);
    APPLY_DBESSEL(g_Bt, bdot_index); // bessel filter
    // APPLY_FBESSEL(g_Bt, bdot_index); // bessel filter
    // printf("readSensors: m0: %d m1: %d Btx: %f Bty: %f Btz: %f\n", m0, m1, x_g_Bt[bdot_index], y_g_Bt[bdot_index], z_g_Bt[bdot_index]);
    getOmega();
    getSVec();
    // log data
    // check if any of the values are NaN. If so, return -1
    // the NaN may stem from Bdot = 0, which may stem from the fact that during sunpointing
    // B may align itself with Z/ω
    if (isnan(x_g_B[mag_index]))
        return -1;
    if (isnan(y_g_B[mag_index]))
        return -1;
    if (isnan(z_g_B[mag_index]))
        return -1;

    if (isnan(x_g_W[omega_index]))
        return -1;
    if (isnan(y_g_W[omega_index]))
        return -1;
    if (isnan(z_g_W[omega_index]))
        return -1;

    if (isnan(x_g_S[sol_index]))
        return -1;
    if (isnan(y_g_S[sol_index]))
        return -1;
    if (isnan(z_g_S[sol_index]))
        return -1;
    return status;
}
#endif