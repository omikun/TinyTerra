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
                grav.velocity += NextVelocity(dt, grav);
                //next position MUST be overwritten
                //trajectory calculation reads and clobbers this value
                grav.nextPosition = grav.transform.position + dt * grav.velocity;
            }
            foreach (var grav in gravities)
            {
                grav.transform.position = grav.nextPosition;
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

    public void UpdateTrajectories() 
    {
        float dt = 0.2f;
        //init
        foreach (var grav in gravities)
        {
            grav.nextVelocity = grav.velocity;
            grav.simPosition = grav.transform.position;
            grav.nextPosition = grav.transform.position;
            grav.nextPositions.Clear();
            grav.nextPositions.Add(grav.transform.position);
        }

        for (int i = 0; i < updateInterval * 500; i++)
        {
            //calculate velocity and position for next iteration
            foreach (var grav in gravities)
            {
                grav.nextVelocity += SimulateNextVelocity(dt, grav);
                grav.nextPosition += dt * grav.nextVelocity;
            }
            //update next to current
            foreach (var grav in gravities)
            {
                grav.nextPositions.Add(grav.nextPosition - grav.soiGrav.nextPosition + grav.soiGrav.transform.position);
                grav.simPosition = grav.nextPosition;
            }
        }

        //update trajectory graphics
        foreach (var grav in gravities)
        {
            grav.UpdateTrajectory(grav.nextPositions);
        }
    }
    Vector3 SimulateNextVelocity(float dt, gravity bodyA) 
    {
        Vector3 velocity = new Vector3(0, 0, 0);
        //pull ship toward gravity
        foreach (var grav in gravities) 
        {
            if (grav.gameObject.GetInstanceID() == bodyA.gameObject.GetInstanceID()) 
            { 
                continue;
            }
            
            var towardsGrav = grav.simPosition - bodyA.simPosition;
            var normalizedTowardGrav = towardsGrav.normalized;
            var sqradius = towardsGrav.sqrMagnitude;
            velocity += normalizedTowardGrav * g * grav.mass * dt * (1 / sqradius);
        }
        return velocity;
    }

    Vector3 NextVelocity(float dt, gravity bodyA) 
    {
        Vector3 velocity = new Vector3(0, 0, 0);
        //pull ship toward gravity
        foreach (var grav in gravities) 
        {
            if (grav.gameObject.GetInstanceID() == bodyA.gameObject.GetInstanceID()) 
            { 
                continue;
            }

            var towardsGrav = grav.transform.position - bodyA.transform.position;
            var normalizedTowardGrav = towardsGrav.normalized;
            var sqradius = towardsGrav.sqrMagnitude;
            velocity += normalizedTowardGrav * g * grav.mass * dt * (1 / sqradius);
        }
        return velocity;
    }
}
