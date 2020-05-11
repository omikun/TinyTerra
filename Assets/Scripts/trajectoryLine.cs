using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class trajectoryLine : MonoBehaviour
{
    LineRenderer lines;
    public Color c1 = Color.yellow;
    public Color c2 = Color.red;
    // Start is called before the first frame update
    void Start()
    {
        gameObject.AddComponent<LineRenderer>();
        lines = gameObject.GetComponent<LineRenderer>();
        lines.widthMultiplier = .2f;
        float alpha = 1.0f;
        Gradient gradient = new Gradient();
        gradient.SetKeys(
            new GradientColorKey[] { new GradientColorKey(c1, 0.0f), new GradientColorKey(c2, 1.0f) },
            new GradientAlphaKey[] { new GradientAlphaKey(alpha, 0.0f), new GradientAlphaKey(alpha, 1.0f) }
        );
        lines.colorGradient = gradient;
    }

    public void UpdateTrajectory(List<Vector3> points)
    {
        lines.positionCount = points.Count;
        lines.SetPositions(points.ToArray());
    }
    // Update is called once per frame
    void Update()
    {
        
    }
}
