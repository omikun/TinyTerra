using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TimeManager : MonoBehaviour
{
    private float fixedDeltaTime;
    [Range(.5f, 3.0f)]
    public float timeScale = 1.0f;
    // Start is called before the first frame update
    void Start()
    {
        this.fixedDeltaTime = Time.fixedDeltaTime;
    }

    // Update is called once per frame
    void Update()
    {
        Time.timeScale = timeScale*timeScale;
        Time.fixedDeltaTime = this.fixedDeltaTime * Time.timeScale;
    }
}
