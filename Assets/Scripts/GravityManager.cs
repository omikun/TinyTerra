using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GravityManager : MonoBehaviour
{
    // Start is called before the first frame update
    gravity[] gravities;
    public float g;
    public bool pause = false;

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
            print (grav.gameObject.name + "->" + grav.soiGrav.gameObject.name + " old v: " + oldV.ToString() + " new v: " + grav.rb.velocity.ToString());
            //print ("toSOI: " + toSOI.ToString() + " velVector: " + velVector.ToString() + "velV.norm: " + velVector.normalized.ToString() + " speed: " + speed.ToString());
        }
        //find chain of soigravs
        
        foreach (var grav in gravities)
        {
            if (IsSame(grav, grav.soiGrav)) { continue; }
            grav.rb.velocity += grav.soiGrav.rb.velocity;
            print(grav.gameObject.name + " v: " + grav.rb.velocity.ToString());
        }
    }

    bool IsSame(gravity a, gravity b )
    {
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

    int updateInterval = 1;
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
            UpdateTrajectories();
        }
        frameCounter += 1;
    }

    public bool useRK4 = false;
    public void UpdateTrajectories() 
    {
        float dt = 0.2f;
        int interval = 100;
        if (useRK4) {
            dt = .6f;
            interval = 50;
        }
        //init
        foreach (var grav in gravities)
        {
            grav.nextVelocity = grav.rb.velocity;
            grav.simPosition = grav.transform.position;
            grav.nextPosition = grav.transform.position;
            grav.nextPositions.Clear();
            grav.nextPositions.Add(grav.transform.position);
            grav.runSimulation = true;
        }

        for (int i = 0; i < updateInterval * interval; i++)
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
                //TODO
                // print if SOI is not grav.soiGrav
                //future pos - future soi pos + current soi pos
                WillCollideWithAny(grav);
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
            
            var towardsGrav = grav.simPosition - bodyA.simPosition;
            var normalizedTowardGrav = towardsGrav.normalized;
            var sqradius = towardsGrav.sqrMagnitude;
            if (useRK4) 
            {
                var gm = g * grav.rb.mass;
                var v = OrbitUtil.rk4(dt, bodyA.simPosition, bodyA.nextVelocity, grav.simPosition, gm, Vector3.zero);
                velocity += v - bodyA.nextVelocity;
            } else {
                var accel = normalizedTowardGrav * g * grav.rb.mass / sqradius;
                velocity += accel * dt;
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

            var towardsGrav = grav.transform.position - bodyA.transform.position;
            var normalizedTowardGrav = towardsGrav.normalized;
            var sqradius = towardsGrav.sqrMagnitude;
            var accel = normalizedTowardGrav * g * grav.rb.mass / sqradius;
            velocity += accel * dt;
        }
        return velocity;
    }
}
