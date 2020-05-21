//from https://github.com/omikun/OrbitalVector/blob/master/Assets/Scripts/Orbit/OrbitalElements.cs
using UnityEngine;
using System.Collections;

using System;
using System.Collections.Generic;

public class OrbitalElements
{
    public float sma;
    public float ecc;
    public float inc;
    public float lan;
    public float aop;
    public float tra;
    public float computeTime;
    public void print()
    {
        Debug.Log("sma: " + OVTools.FormatDistance((float)sma));
        Debug.Log("ecc: " + ecc);
        Debug.Log("inc: " + inc);
        Debug.Log("lan: " + lan);
        Debug.Log("aop: " + aop);
        Debug.Log("tra: " + tra);
        Debug.Log("compTime: " + OVTools.FormatDistance((float)computeTime));
    }

    public float GetPeriod(float parentGM)
    {
        return 2 * Mathf.PI * Mathf.Sqrt(Mathf.Pow(sma, 3) / parentGM);
    }
    public void SetComputeTime(float t)
    {
        computeTime = t;
    }
    public OrbitalElements CopyOE()
    {
        OrbitalElements ret = new OrbitalElements();
        ret.sma = sma;
        ret.lan = lan;
        ret.inc = inc;
        ret.ecc = ecc;
        ret.tra = tra;
        ret.aop = aop;
        ret.computeTime = computeTime;
        return ret;
    }
    /*
    public Vector3 GetRAtTime(float absTime, float parentGM)
    {
        var tempOE = CopyOE();
        tempOE.tra = Program.anomalyAfterTime(parentGM, tempOE, absTime - computeTime);
        return Util.oe2r(parentGM, tempOE);
    }

    public Vector3 GetRdAtTime(float absTime, float parentGM)
    {
        var tempOE = CopyOE();
        tempOE.tra = Program.anomalyAfterTime(parentGM, tempOE, absTime - computeTime);
        return Util.oe2rd(parentGM, tempOE);
    }
    // */


}
public class Program
{

    /**
* Computes the delta-time until a particular true-anomaly is reached.
* @param       grav_param      The gravitational parameter of the two-body
system.
* @param       oe      The orbital elements of the body of interest.
* @param       true_anomaly    The true anomaly at the final time.
* @returns     The amount of time until the true anomaly has been reached.
*/
    public static double Atanh(double x)
    {
        return (Math.Log(1 + x) - Math.Log(1 - x)) / 2;
    }
    public static float timeUntilAnomaly(float grav_param, OrbitalElements oe, float
    true_anomaly)
    {
        var ret = _timeUntilAnomaly(grav_param, oe, (double)true_anomaly);
        return (float)ret;
    }
    public static double _timeUntilAnomaly(float grav_param, OrbitalElements oe, double
    true_anomaly)
    {
        var sma = (double)oe.sma;
        var tra = (double)oe.tra;
        var ecc = (double)oe.ecc;

        if (sma > 0.0)
        {
            double ecc_anomaly_i =
2.0 * Math.Atan(Math.Sqrt((1.0 - ecc) / (1.0 + ecc)) * Math.Tan(tra / 2.0));
            double ecc_anomaly_f =
2.0 * Math.Atan(Math.Sqrt((1 - ecc) / (1 + ecc)) * Math.Tan(true_anomaly / 2));

            if (ecc_anomaly_f < ecc_anomaly_i) ecc_anomaly_f +=
2.0 * Math.PI;

            double mean_anomaly_i =
ecc_anomaly_i - ecc * Math.Sin(ecc_anomaly_i);
            double mean_anomaly_f =
ecc_anomaly_f - ecc * Math.Sin(ecc_anomaly_f);

            return
Math.Sqrt(Math.Pow(sma, 3.0) / grav_param) * (mean_anomaly_f - mean_anomaly_i);
        }
        else
        {
            double hyp_anomaly_i = 2.0 * Atanh(Math.Sqrt((ecc - 1) / (1 + ecc)) * Math.Tan(tra / 2));
            double hyp_anomaly_f = 2.0 * Atanh(Math.Sqrt((ecc - 1) / (1 + ecc)) * Math.Tan(true_anomaly / 2));

            double mean_anomaly_i =
-hyp_anomaly_i + ecc * Math.Sinh(hyp_anomaly_i);
            double mean_anomaly_f =
-hyp_anomaly_f + ecc * Math.Sinh(hyp_anomaly_f);

            return
Math.Sqrt(Math.Pow(-sma, 3.0) / grav_param) * (mean_anomaly_f - mean_anomaly_i);
        }
    }


    /**
* Computes the true anomaly after a delta-time has passed.
* @param       grav_param      The gravitational parameter of the two-body
system.
* @param       oe      The orbital elements of the body of interest.
* @param       delta_time      The elapsed time.
* @returns     The true anomaly at the final time.
*/
    /*
            public static float anomalyAfterTime(float grav_param, OrbitalElements oe, float
            delta_time)
            {
                if (oe.sma > 0)
                {
                    float ecc_anomaly_i =
        2 * Mathf.Atan(Mathf.Sqrt((1 - oe.ecc) / (1 + oe.ecc)) * Mathf.Tan(oe.tra / 2));
                    float mean_anomaly_f =
        ecc_anomaly_i - oe.ecc * Mathf.Sin(ecc_anomaly_i) + Mathf.Sqrt(grav_param / Mathf.Pow(oe.sma, 3)) * delta_time;

                    // perform Newton-Raphson iteration to determine the final eccentric anomaly
                    float ecc_anomaly_f = mean_anomaly_f;
                    float error =
        mean_anomaly_f - ecc_anomaly_f + oe.ecc * Mathf.Sin(ecc_anomaly_f);
                    while (Mathf.Abs(error) > 1E-10)
                    {
                        ecc_anomaly_f =
        ecc_anomaly_f - error / (oe.ecc * Mathf.Cos(ecc_anomaly_f) - 1);
                        error =
        mean_anomaly_f - ecc_anomaly_f + oe.ecc * Mathf.Sin(ecc_anomaly_f);
                    }

                    return
        2 * Mathf.Atan(Mathf.Sqrt((1 + oe.ecc) / (1 - oe.ecc)) * Mathf.Tan(ecc_anomaly_f / 2));
                }
                else
                {
                    float hyp_anomaly_i =
        2 * Atanh(Mathf.Sqrt((oe.ecc - 1) / (1 + oe.ecc)) * Mathf.Tan(oe.tra / 2));
                    float mean_anomaly_f =
        -hyp_anomaly_i + oe.ecc * Mathf.Sinh(hyp_anomaly_i) + Mathf.Sqrt(grav_param / Mathf.Pow(-oe.sma, 3)) * delta_time;

                    // perform Newton-Raphson iteration to determine the final eccentric anomaly
                    float hyp_anomaly_f = mean_anomaly_f;
                    float error =
        mean_anomaly_f + hyp_anomaly_f - oe.ecc * Mathf.Sinh(hyp_anomaly_f);
                    while (Mathf.Abs(error) > 1E-10)
                    {
                        hyp_anomaly_f =
        hyp_anomaly_f - error / (oe.ecc * Mathf.Cosh(hyp_anomaly_f) - 1);
                        error =
        mean_anomaly_f + hyp_anomaly_f - oe.ecc * Mathf.Sinh(hyp_anomaly_f);
                    }

                    return
        2 * Mathf.Atan(Mathf.Sqrt((1 + oe.ecc) / (oe.ecc - 1)) * Mathf.Tanh(hyp_anomaly_f / 2));
                }
            }

        // */

}
public class Util
{
    /**
* Converts orbital elements to position and velocity.
* @param grav_param The gravitational parameter of the two-body system.
* @param oe The orbital elements to be converted from.
* @returns A vector corresponding to the concatenated position and
velocity.
*/
    public static Vector3 oe2r(float grav_param, OrbitalElements oe)
    {
        var v = oe2rd(grav_param, oe);
        Vector3 pos = new Vector3((float)v.x, (float)v.y, (float)v.z);
        return pos;
    }
    public static Vector3 oe2rd(float grav_param, OrbitalElements oe)
    {
        // rotation matrix
        float R11 =
Mathf.Cos(oe.aop) * Mathf.Cos(oe.lan) - Mathf.Cos(oe.inc) * Mathf.Sin(oe.aop) * Mathf.Sin(oe.lan);
        float R12 =
-Mathf.Cos(oe.lan) * Mathf.Sin(oe.aop) - Mathf.Cos(oe.inc) * Mathf.Cos(oe.aop) * Mathf.Sin(oe.lan);
        //float R13 = Mathf.Sin(oe.inc)*Mathf.Sin(oe.lan);
        float R21 =
Mathf.Cos(oe.inc) * Mathf.Cos(oe.lan) * Mathf.Sin(oe.aop) + Mathf.Cos(oe.aop) * Mathf.Sin(oe.lan);
        float R22 =
Mathf.Cos(oe.inc) * Mathf.Cos(oe.aop) * Mathf.Cos(oe.lan) - Mathf.Sin(oe.aop) * Mathf.Sin(oe.lan);
        //float R23 = -Mathf.Cos(oe.lan)*Mathf.Sin(oe.inc);
        float R31 = Mathf.Sin(oe.inc) * Mathf.Sin(oe.aop);
        float R32 = Mathf.Cos(oe.aop) * Mathf.Sin(oe.inc);
        //float R33 = Mathf.Cos(oe.inc);

        // semi-latus rectum
        float p = oe.sma * (1 - oe.ecc * oe.ecc);

        // position in the perifocal frame
        //std::vector<float> r_pf(3);
        float[] r_pf = new float[3];
        float r_norm = p / (1 + oe.ecc * Mathf.Cos(oe.tra));
        r_pf[0] = r_norm * Mathf.Cos(oe.tra);
        r_pf[1] = r_norm * Mathf.Sin(oe.tra);
        r_pf[2] = 0;

        // velocity in the perifocal frame
        //std::vector<float> v_pf(3);
        float[] v_pf = new float[3];
        v_pf[0] = Mathf.Sqrt(grav_param / p) * -Mathf.Sin(oe.tra);
        v_pf[1] = Mathf.Sqrt(grav_param / p) * (oe.ecc + Mathf.Cos(oe.tra));
        v_pf[2] = 0;

        // rotate the position and velocity into the body-fixed inertial frame
        //std::vector<float> rv(6);
        VectorD r = new VectorD();
        r.Resize(3);
        r[0] = R11 * r_pf[0] + R12 * r_pf[1] /*+R13*r_pf[2]*/;
        r[1] = R21 * r_pf[0] + R22 * r_pf[1] /*+R23*r_pf[2]*/;
        r[2] = R31 * r_pf[0] + R32 * r_pf[1] /*+R33*r_pf[2]*/;

        Vector3 pos = new Vector3(r[0], r[1], r[2]);
        return pos;
    }
    public static void oe2rv(float grav_param, OrbitalElements oe, out Vector3 r, out Vector3 v)
    {
        var rv = oe2rv(grav_param, oe);
        r = new Vector3(rv[0], rv[1], rv[2]);
        v = new Vector3(rv[3], rv[4], rv[5]);
    }
    public static VectorD oe2rv(float grav_param, OrbitalElements oe)
    {
        // rotation matrix
        float R11 =
Mathf.Cos(oe.aop) * Mathf.Cos(oe.lan) - Mathf.Cos(oe.inc) * Mathf.Sin(oe.aop) * Mathf.Sin(oe.lan);
        float R12 =
-Mathf.Cos(oe.lan) * Mathf.Sin(oe.aop) - Mathf.Cos(oe.inc) * Mathf.Cos(oe.aop) * Mathf.Sin(oe.lan);
        //float R13 = Mathf.Sin(oe.inc)*Mathf.Sin(oe.lan);
        float R21 =
Mathf.Cos(oe.inc) * Mathf.Cos(oe.lan) * Mathf.Sin(oe.aop) + Mathf.Cos(oe.aop) * Mathf.Sin(oe.lan);
        float R22 =
Mathf.Cos(oe.inc) * Mathf.Cos(oe.aop) * Mathf.Cos(oe.lan) - Mathf.Sin(oe.aop) * Mathf.Sin(oe.lan);
        //float R23 = -Mathf.Cos(oe.lan)*Mathf.Sin(oe.inc);
        float R31 = Mathf.Sin(oe.inc) * Mathf.Sin(oe.aop);
        float R32 = Mathf.Cos(oe.aop) * Mathf.Sin(oe.inc);
        //float R33 = Mathf.Cos(oe.inc);

        // semi-latus rectum
        float p = oe.sma * (1 - oe.ecc * oe.ecc);

        // position in the perifocal frame
        //std::vector<float> r_pf(3);
        float[] r_pf = new float[3];
        float r_norm = p / (1 + oe.ecc * Mathf.Cos(oe.tra));
        r_pf[0] = r_norm * Mathf.Cos(oe.tra);
        r_pf[1] = r_norm * Mathf.Sin(oe.tra);
        r_pf[2] = 0;

        // velocity in the perifocal frame
        //std::vector<float> v_pf(3);
        float[] v_pf = new float[3];
        v_pf[0] = Mathf.Sqrt(grav_param / p) * -Mathf.Sin(oe.tra);
        v_pf[1] = Mathf.Sqrt(grav_param / p) * (oe.ecc + Mathf.Cos(oe.tra));
        v_pf[2] = 0;

        // rotate the position and velocity into the body-fixed inertial frame
        //std::vector<float> rv(6);
        VectorD rv = new VectorD();
        rv.Resize(6);
        rv[0] = R11 * r_pf[0] + R12 * r_pf[1] /*+R13*r_pf[2]*/;
        rv[1] = R21 * r_pf[0] + R22 * r_pf[1] /*+R23*r_pf[2]*/;
        rv[2] = R31 * r_pf[0] + R32 * r_pf[1] /*+R33*r_pf[2]*/;
        rv[3] = R11 * v_pf[0] + R12 * v_pf[1] /*+R13*v_pf[2]*/;
        rv[4] = R21 * v_pf[0] + R22 * v_pf[1] /*+R23*v_pf[2]*/;
        rv[5] = R31 * v_pf[0] + R32 * v_pf[1] /*+R33*v_pf[2]*/;

        return rv;
    }


    /**
     * Converts position and velocity to orbital elements
     * @param grav_param The gravitational parameter of the two-body system.
     * @param rv The concatenated position and velocity to convert from.
     * @returns The orbital elements.
     */
    public static OrbitalElements rv2oe(float grav_param, VectorD rv)
    {

        OrbitalElements oe = new OrbitalElements();

        // Semi-major Axis : Vis-viva Equation
        oe.sma = 1.0f / (
        2.0f / Mathf.Sqrt(rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2])
        - (rv[3] * rv[3] + rv[4] * rv[4] + rv[5] * rv[5]) / grav_param
        );

        // Angular Momentum
        //List<float> h = new List<float>(3);
        float[] h = new float[3];
        h[0] = rv[1] * rv[5] - rv[2] * rv[4];
        h[1] = rv[2] * rv[3] - rv[0] * rv[5];
        h[2] = rv[0] * rv[4] - rv[1] * rv[3];

        // Norm of position
        float r = Mathf.Sqrt(rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]);

        // Eccentricity Vector :  e = v x h / mu - r/|r|
        //std::vector<float> e(3);
        float[] e = new float[3];
        e[0] = (rv[4] * h[2] - rv[5] * h[1]) / grav_param - rv[0] / r;
        e[1] = (rv[5] * h[0] - rv[3] * h[2]) / grav_param - rv[1] / r;
        e[2] = (rv[3] * h[1] - rv[4] * h[0]) / grav_param - rv[2] / r;

        // Eccentricity
        oe.ecc = Mathf.Sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);

        // Inclination
        oe.inc = Mathf.Acos(h[2] / Mathf.Sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]));

        // Ascending Node Direction (In x-y plane)
        //std::vector<float> n(2);
        float[] n = new float[2];
        n[0] = -h[1] / Mathf.Sqrt(h[0] * h[0] + h[1] * h[1]);
        n[1] = h[0] / Mathf.Sqrt(h[0] * h[0] + h[1] * h[1]);
        float n_norm = Mathf.Sqrt(n[0] * n[0] + n[1] * n[1]);

        // Longitude of the Ascending Node
        oe.lan = Mathf.Acos(n[0]) / n_norm;
        if (n[1] < 0.0) oe.lan = 2 * Mathf.PI - oe.lan;

        // Argument of Periapsis
        oe.aop = Mathf.Acos((n[0] * e[0] + n[1] * e[1]) / (n_norm * oe.ecc));
        if (e[2] < 0.0) oe.aop = 2 * Mathf.PI - oe.aop;

        // True Anomaly
        oe.tra = Mathf.Acos((rv[0] * e[0] + rv[1] * e[1] + rv[2] * e[2]) / (r * oe.ecc));
        if (rv[0] * rv[3] + rv[1] * rv[4] + rv[2] * rv[5] < 0.0) oe.tra = 2 * Mathf.PI - oe.tra;

        return oe;
    }
    public static VectorD convertToRv(ref Vector3 pos, ref Vector3 vel)
    {
        VectorD params_ = new VectorD();
        params_.Resize(6);
        Debug.Assert(params_.Count == 6);
        //2013 Orbital Vector was written with z being up
        //unity y is up, so must swap axis
        params_[0] = pos.x;
        params_[1] = pos.y;
        params_[2] = pos.z;
        params_[3] = vel.x;
        params_[4] = vel.y;
        params_[5] = vel.z;
        return params_;
    }
    public static VectorD convertToRv(ref float[] pos, ref float[] vel)
    {
        Debug.Assert(vel.Length == 3);
        Debug.Assert(pos.Length == 3);

        VectorD params_ = new VectorD();
        params_.Resize(6);
        Debug.Assert(params_.Count == 6);
        //2013 Orbital Vector was written with z being up
        //unity y is up, so must swap axis
        params_[0] = pos[0];
        params_[1] = pos[1];
        params_[2] = pos[2];
        params_[3] = vel[0];
        params_[4] = vel[1];
        params_[5] = vel[2];
        return params_;
    }
    public static VectorD convertToParams(float[] parentPos, float gm, float[] accel)
    {
        VectorD params_ = new VectorD();
        params_.Resize(7);
        params_[0] = parentPos[0];
        params_[1] = parentPos[1];
        params_[2] = parentPos[2];
        params_[3] = gm;
        params_[4] = accel[0];
        params_[5] = accel[1];
        params_[6] = accel[2];
        return params_;
    }

    public static VectorD convertToParams(Vector3 parentPos, float gm, Vector3 accel)
    {
        VectorD params_ = new VectorD();
        params_.Resize(7);
        //parent pos, gm, accel
        params_[0] = parentPos.x;
        params_[1] = parentPos.y;
        params_[2] = parentPos.z;
        params_[3] = gm;
        params_[4] = accel.x;
        params_[5] = accel.y;
        params_[6] = accel.z;
        return params_;
    }
    
    public static Vector3 rk4(float dt, Vector3 pos, Vector3 vel, Vector3 parentPos, float gm, Vector3 accel)
    {
        var t0 = (float)dt;
        var x0 = convertToRv(ref pos, ref vel);
        var params_ = convertToParams(parentPos, gm, accel);

        var result = rungeKutta4(0, t0, x0, params_);
        var ret = new Vector3((float)result[3], (float)result[4], (float)result[5]);
        return ret;
    }
    public static VectorD forwardEuler(float t0, float dt, VectorD x0, VectorD params_)
    {
        return x0 + dt * dynamics(t0, x0, params_);
    }

    public static VectorD rungeKutta4(float t0, float dt, VectorD x0, VectorD params_)
    {
        VectorD k1 = dynamics(t0, x0, params_);
        VectorD k2 = dynamics(t0 + dt / 2, x0 + k1 * dt / 2, params_);
        VectorD k3 = dynamics(t0 + dt / 2, x0 + k2 * dt / 2, params_);
        VectorD k4 = dynamics(t0 + dt, x0 + k3 * dt, params_);
        return x0 + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    }

    public static VectorD dynamics(float t, VectorD x, VectorD params_)
    {
        Debug.Assert(x.Count == 6);
        Debug.Assert(params_.Count == 7);

        VectorD dx = new VectorD();
        dx.Resize(6);

        //the derivative of position is velocity
        dx[0] = x[3];
        dx[1] = x[4];
        dx[2] = x[5];

        //distance between the two bodies
        float distance = 0;
        distance += (x[0] - params_[0]) * (x[0] - params_[0]);
        distance += (x[1] - params_[1]) * (x[1] - params_[1]);
        distance += (x[2] - params_[2]) * (x[2] - params_[2]);
        distance = Mathf.Sqrt(distance);

        //acceleration due to gravity
        dx[3] = -params_[3] * (x[0] - params_[0]) / Mathf.Pow(distance, 3);
        dx[4] = -params_[3] * (x[1] - params_[1]) / Mathf.Pow(distance, 3);
        dx[5] = -params_[3] * (x[2] - params_[2]) / Mathf.Pow(distance, 3);

        //perturbing acceleration
        dx[3] += params_[4];
        dx[4] += params_[5];
        dx[5] += params_[6];

        return dx;

    }
}