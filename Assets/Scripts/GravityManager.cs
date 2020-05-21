using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GravityManager : MonoBehaviour
{
    // Start is called before the first frame update
    gravity[] gravities;
    public float g;
    public bool pause = false;
    float maxProjectionTime = 300;
    float standardMaxProjectionTime = 300;

    //record all gravities
    void Start()
    {
        gravities = GameObject.FindObjectsOfType<gravity>();

        
        //foreach (var grav in gravities)
        //{
        //    SetSOIFor(grav);
        //}

        //set all grav into circular orbits around SOI
        foreach (var grav in gravities)
        {
            if (IsSame(grav, grav.soiGrav)) { continue; }
            //get vector to soi
            var toSOI = grav.transform.position - grav.soiGrav.transform.position;
            //cross with up to get velocity vector
            var velVector = Vector3.Cross(toSOI, Vector3.up);
            //get speed
            var speed = Mathf.Pow(g * grav.soiGrav.rb.mass / toSOI.magnitude, 0.5f);
            var oldV = grav.rb.velocity;
            var v = velVector.normalized * speed;
            grav.rb.velocity = v;
            print (grav.gameObject.name + " old v: " + oldV.ToString());
            print (grav.gameObject.name + "->" + grav.soiGrav.gameObject.name + " new v: " + grav.rb.velocity.ToString());
            //print ("toSOI: " + toSOI.ToString() + " velVector: " + velVector.ToString() + "velV.norm: " + velVector.normalized.ToString() + " speed: " + speed.ToString());
        }
        //find chain of soigravs
        
        foreach (var grav in gravities)
        {
            if (IsSame(grav, grav.soiGrav)) { continue; }
            grav.nextVelocity = Vector3.zero;
            var tempGrav = grav;
            print ("Adding for " + tempGrav.gameObject.name);
            while (!IsSame(tempGrav, tempGrav.soiGrav))
            {
                grav.nextVelocity += tempGrav.soiGrav.rb.velocity;
                //print("\t+: " + tempGrav.soiGrav.gameObject.name + " v: " + tempGrav.soiGrav.nextVelocity.ToString());
                tempGrav = tempGrav.soiGrav;
            }
            //grav.rb.velocity += grav.soiGrav.rb.velocity;
            //print(grav.gameObject.name + " v: " + grav.rb.velocity.ToString());
        }
        
        foreach (var grav in gravities)
        {
            grav.rb.velocity += grav.nextVelocity;
            print (grav.gameObject.name + "->" + grav.soiGrav.gameObject.name + " new v: " + grav.rb.velocity.ToString());
        }
    }

    bool IsSame(gravity a, gravity b )
    {
        if (a == null || b == null) { return false; }
        return (a.gameObject.GetInstanceID() == b.gameObject.GetInstanceID());
    }

    struct SOI
    {
        public float dist;
        public gravity grav;
    }
    void UpdateSOI(gravity grav)
    {
        var soi = GetSOI(grav);

        if (soi.grav == null)
        {
            soi.grav = grav.soiGrav;
        }
        else
        {
            if (!IsSame(soi.grav, grav.soiGrav))
            {
                print(grav.gameObject.name + "'s old soi: " + grav.soiGrav.gameObject.name + " new soi: " + soi.grav.gameObject.name);
            }
            grav.soiGrav = soi.grav;
        }
    }
    SOI GetSOI(gravity thisGrav)
    {
        var ret = new SOI{ dist = float.MaxValue, grav = null};
        foreach (var grav in gravities)
        {
            //find nearest grav
            if (IsSame(grav, thisGrav)) { continue; }
            if (grav.rb.mass < thisGrav.rb.mass) { continue; }

            var soiDist = GetSOIDist(grav, thisGrav);
            if (soiDist < ret.dist)
            {
                ret.dist = soiDist;
                ret.grav = grav;
            }
        }
        return ret;
    }

    gravity GetNearestSOIGrav(gravity notThis, Vector3 position)
    {
        var largestAccel = float.MinValue;
        gravity ret = null;
        foreach (var a in gravities)
        {
            if (IsSame(notThis, a)) { continue; }
            if (a.rb.mass < notThis.rb.mass) { continue; }
            var accel = g * a.rb.mass / (position - a.nextPosition).sqrMagnitude;
            if (largestAccel < accel)
            {
                largestAccel = accel;
                ret = a;
            }
        }
        return ret;
    }
    void SetSOIFor(gravity thisGrav)
    {
        var soi = GetSOI(thisGrav);
        if (soi.grav == null)
        {
            print("error! couldn't find SOI for this grav");
            return;
        }
        thisGrav.soiGrav = soi.grav;
    }

    float GetSOIDist(gravity a, gravity b) {
        var m1 = a.rb.mass/b.rb.mass;
        var m2 = 1/m1;
        var minM = Mathf.Min(m1, m2);
        var dist = (a.transform.position - b.transform.position).magnitude; 
        return dist * Mathf.Pow(minM, 2/5f);
    }

    int updateInterval = 80;
    int frameCounter = 0;
    // Update is called once per frame
    void Update()
    {
        var dt = Time.deltaTime;
        //compute next positions for all gravities
        if (!pause)
        {
            foreach (var grav in gravities)
            {
                grav.rb.velocity += NextVelocity(dt, grav);
                //UpdateSOI(grav);
            }
        }

        //update positions to next positions for all gravities
        if (frameCounter % updateInterval == 0) 
        {
            frameCounter = 0;
            foreach (var grav in gravities)
            {
                print(grav.gameObject.name + "->" + grav.soiGrav.gameObject.name + " new v: " + grav.rb.velocity.magnitude.ToString());
            }
        }
            UpdateTrajectories();
        frameCounter += 1;
    }

    [Range(.1f, 2f)]
    public float updateDt = .2f;
    [Range(10, 1000)]
    public int numUpdates = 100;
    public bool useRK4 = false;
    [Range(.1f, 10f)]
    public float rk4UpdateDt = .2f;
    [Range(50, 1000)]
    public int rk4NumUpdates = 100;
    public void UpdateTrajectories() 
    {
        float dt = updateDt;
        int interval = numUpdates;
        if (useRK4) {
            dt = rk4UpdateDt;
            interval = rk4NumUpdates;
        }
        maxProjectionTime = dt * interval;
        standardMaxProjectionTime = dt * interval;
        //init
        foreach (var grav in gravities)
        {
            grav.nextVelocity = grav.rb.velocity;
            grav.simPosition = grav.transform.position;
            grav.nextPosition = grav.transform.position;
            grav.nextPositions.Clear();
            grav.nextPositions.Add(grav.transform.position);
            grav.runSimulation = true;
            grav.simSoiGrav = grav.soiGrav;
            grav.simPoint = 0;
            grav.maxProjectionTime = standardMaxProjectionTime;
        }

        for (int i = 0; i < maxProjectionTime; i++)
        {
            //calculate velocity and position for next iteration
            foreach (var grav in gravities)
            {
                if (grav.runSimulation == false) { continue; }
                grav.nextVelocity += SimulateNextVelocity(dt, grav);
                grav.nextPosition += dt * grav.nextVelocity;
            }
            //update next to current
            foreach (var grav in gravities)
            {
                var soi = GetNearestSOIGrav(grav, grav.nextPosition);
                if (soi && IsSame(soi, grav.simSoiGrav) == false)
                {
                    grav.simSoiGrav = soi;
                    grav.simPoint = dt * i;
                    LimitToPeriod(grav, i * dt);
                }
                if (soi == null)
                {
                    soi = grav.soiGrav;
                }
                //TODO
                // print if SOI is not grav.soiGrav
                //future pos - future soi pos + current soi pos
                WillCollideWithAny(grav);
                ShutDown(grav, i * dt);
                if (grav.runSimulation == false) { continue; }

                grav.nextPositions.Add(grav.nextPosition - soi.nextPosition + soi.transform.position);
                grav.simPosition = grav.nextPosition;
            }
        }

        //update trajectory graphics
        foreach (var grav in gravities)
        {
            grav.UpdateTrajectory(grav.nextPositions);
        }
    }

    void ShutDown(gravity grav, float t)
    {
        if (grav.maxProjectionTime <= t) 
        {
            grav.runSimulation = false;
        }
    }
    //if period is too long, stop displaying this thing
    void LimitToPeriod(gravity grav, float currentTime)
    {
        var period = GetPeriod(grav, grav.nextPosition, grav.nextVelocity);

        if (currentTime + period > maxProjectionTime && period < standardMaxProjectionTime)
        {
            maxProjectionTime = currentTime + period;
            grav.maxProjectionTime = maxProjectionTime;
        }
    }
    float GetPeriod(gravity grav, Vector3 pos, Vector3 vel)
    {
        var rv = Util.convertToRv(ref pos, ref vel);
        
        var gm = g * grav.soiGrav.rb.mass;
        var oe = Util.rv2oe(gm, rv);

        //find period of gravity around gravwell
        if (oe.sma < .1f)
        {
            return 0;
        } else {
            var p = oe.GetPeriod(gm);
            if (p > maxProjectionTime)
            {
                return 0;
            }
            return p;
        }
    }
    bool WillCollideWithAny(gravity a) {
        foreach (var grav in gravities)
        {
            if (IsSame(a, grav)) { continue; }
            var boom = IsCollidingWith(grav.nextPosition, a.nextPosition, grav.radius + a.radius);
            if (boom) 
            {
                if (grav.rb.mass > a.rb.mass)
                {
                    a.runSimulation = false;
                } else {
                    grav.runSimulation = false;
                }
                return true;
            }
        }
        return false;
    }

    bool IsCollidingWith(Vector3 a, Vector3 b, float minDist)
    {
        var dist = (a - b).magnitude;
        return dist <= minDist;
    }
    bool IsCollidingWith(gravity a, gravity b) 
    {
        if (!a || !b) { return false; }

        var dist = (a.transform.position - b.transform.position).magnitude;
        var minDist = a.radius + b.radius;
        return dist <= minDist;
    }
    Vector3 SimulateNextVelocity(float dt, gravity bodyA) 
    {
        Vector3 velocity = new Vector3(0, 0, 0);
        //pull ship toward gravity
        foreach (var grav in gravities) 
        {
            if (IsSame(grav, bodyA)) 
            { 
                continue;
            }
            
            if (useRK4) 
            {
                var gm = g * grav.rb.mass;
                var v = Util.rk4(dt, bodyA.simPosition, bodyA.nextVelocity, grav.simPosition, gm, Vector3.zero);
                velocity += v - bodyA.nextVelocity;
            } else {
                velocity += CalculateVelocity(dt, bodyA.simPosition, grav.simPosition, grav.rb.mass);
            }
        }
        return velocity;
    }

    Vector3 NextVelocity(float dt, gravity bodyA) 
    {
        Vector3 velocity = new Vector3(0, 0, 0);
        //pull ship toward gravity
        foreach (var grav in gravities) 
        {
            if (IsSame(grav, bodyA)) 
            { 
                continue;
            }

            velocity += CalculateVelocity(dt, bodyA.transform.position, grav.transform.position, grav.rb.mass);
        }
        return velocity;
    }

    Vector3 CalculateVelocity(float dt, Vector3 pos, Vector3 parentPos, float parentMass)
    {
        var towardsGrav = parentPos- pos;
        var normalizedTowardGrav = towardsGrav.normalized;
        var sqradius = towardsGrav.sqrMagnitude;
        var accel = normalizedTowardGrav * g * parentMass / sqradius;
        //different gravitational curves:
        //1/(1+e^-(-x+6)), from 0 to 12
        //or
        //{ ((cos(x*pi/5)+1)/2)^2, x<3 }, 1/x^2 otherwise
        if (false && sqradius < 9) 
        {
            var radius = Mathf.Sqrt(sqradius);
            accel = normalizedTowardGrav * g * parentMass 
                * Mathf.Pow( (Mathf.Cos(radius * Mathf.PI / 5) + 1) / 2, 2);
        }
        return accel * dt;
    }
}
