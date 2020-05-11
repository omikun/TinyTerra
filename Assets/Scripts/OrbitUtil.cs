using System.Collections;
using System;
using System.Collections.Generic;
using UnityEngine;

public class VectorD : List<float>
{
    public void Print(string desc)
    {
        return;
        System.Text.StringBuilder sb = new System.Text.StringBuilder();
        sb.Append(desc);
        this.ForEach(x => sb.Append(x + " "));
        Debug.Log(sb.ToString());
    }
    public void Resize(int count)
    {
        count = count - this.Count;
        for (int i = 0; i < count; i++)
        {
            this.Add(0);
        }
    }
    public static VectorD operator +(VectorD a, VectorD b)
    {

        Debug.Assert(a.Count == b.Count);

        VectorD result = new VectorD();

        for (int i = 0; i < a.Count; ++i)
        {
            result.Add(a[i] + b[i]);
        }

        return result;
    }
    public static VectorD operator *(float a, VectorD b)
    {
        VectorD result = new VectorD();

        for (int i = 0; i < b.Count; ++i)
        {
            result.Add(a * b[i]);
        }

        return result;
    }

    public static VectorD operator *(VectorD a, float b)
    {
        return b * a;
    }

    public static VectorD operator /(VectorD a, float b)
    {
        VectorD result = new VectorD();

        for (int i = 0; i < a.Count; ++i)
        {
            result.Add(a[i] / b);
        }

        return result;
    }
}

class OrbitUtil 
{
    public static Vector3 rk4(float dt, Vector3 pos, Vector3 vel, Vector3 parentPos, float gm, Vector3 accel)
    {
        var t0 = (float)dt;
        VectorD x0 = new VectorD();
        x0.Resize(6);
        //pos,vel
        x0[0] = pos.x;
        x0[1] = pos.y;
        x0[2] = pos.z;
        x0[3] = vel.x;
        x0[4] = vel.y;
        x0[5] = vel.z;
        
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

        var result = rungeKutta4(0, t0, x0, params_);
        var ret = new Vector3((float)result[3], (float)result[4], (float)result[5]);
        return ret;
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