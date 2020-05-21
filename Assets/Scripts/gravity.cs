using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class gravity : MonoBehaviour
{
    public Rigidbody rb;
    public float radius;
    public GameObject soi;
    public gravity soiGrav;
    public gravity simSoiGrav; //soi grav at simPosition
    public float simPoint; //time at which simSoiGrav changes
    public float maxProjectionTime;
    public Vector3 velocity;
    public Vector3 nextVelocity = new Vector3(1, 0 , 0);
    public Vector3 nextPosition;
    public Vector3 simPosition;
    public bool runSimulation = false;
    public List<Vector3> nextPositions = new List<Vector3>();
    LineRenderer lines;
    public Color c1 = Color.yellow;
    public Color c2 = Color.red;
    // Start is called before the first frame update
    void Awake()
    {
        AddLineRenderer();
        soiGrav = soi.GetComponent<gravity>();
        rb = GetComponent<Rigidbody>();
        //rb.velocity = velocity;

        var planet = GetComponent<Planet>();
        if (planet) {
            radius = planet.shapeSettings.planetRadius;
        }
    }

    void AddLineRenderer()
    {
        gameObject.AddComponent<LineRenderer>();
        lines = gameObject.GetComponent<LineRenderer>();
        lines.useWorldSpace = true;
        lines.widthMultiplier = 0.03f;
        float alpha = 1.0f;
        Gradient gradient = new Gradient();
        gradient.SetKeys(
            new GradientColorKey[] { new GradientColorKey(c1, 0.0f), new GradientColorKey(c2, 1.0f) },
            new GradientAlphaKey[] { new GradientAlphaKey(alpha, 0.0f), new GradientAlphaKey(alpha, 1.0f) }
        );
        lines.colorGradient = gradient;
        //lines.material.color = c2;
    }

    public bool IsColliding() 
    {
        if (!soiGrav) { return false; }

        var dist = (transform.position - soiGrav.transform.position).magnitude;
        var minDist = radius + soiGrav.radius;
        return dist <= minDist;
    }
    public void UpdateTrajectory(List<Vector3> points)
    {
        lines.positionCount = points.Count;
        lines.SetPositions(points.ToArray());
    }
}
