using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class gravity : MonoBehaviour
{
    public GameObject[] gravities;
    GameObject star;
    public GameObject soi;
    public gravity soiGrav;
    public Vector3 velocity = new Vector3(1, 0 , 0);
    public Vector3 nextVelocity = new Vector3(1, 0 , 0);
    [Range(0.0f, 10.0f)]
    public float velX;
    [Range(0.0f, 10.0f)]
    public float velY;
    public float mass = 10;
    public Vector3 nextPosition;
    public Vector3 simPosition;
    public List<Vector3> nextPositions = new List<Vector3>();
    float g = 100;
    trajectoryLine trajectory;
    LineRenderer lines;
    public Color c1 = Color.yellow;
    public Color c2 = Color.red;
    // Start is called before the first frame update
    void Start()
    {
        AddLineRenderer();
        soiGrav = soi.GetComponent<gravity>();
    }

    [ExecuteInEditMode]
    void AddLineRenderer()
    {
        gameObject.AddComponent<LineRenderer>();
        lines = gameObject.GetComponent<LineRenderer>();
        lines.useWorldSpace = true;
        lines.widthMultiplier = .1f;
        float alpha = 1.0f;
        Gradient gradient = new Gradient();
        gradient.SetKeys(
            new GradientColorKey[] { new GradientColorKey(c1, 0.0f), new GradientColorKey(c2, 1.0f) },
            new GradientAlphaKey[] { new GradientAlphaKey(alpha, 0.0f), new GradientAlphaKey(alpha, 1.0f) }
        );
        lines.colorGradient = gradient;
        //lines.material.color = c2;
    }

    public void UpdateTrajectory(List<Vector3> points)
    {
        lines.positionCount = points.Count;
        lines.SetPositions(points.ToArray());
    }
}