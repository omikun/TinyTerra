using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Indicator : MonoBehaviour
{
    GameObject parent;
    public GameObject target;
    [Range(.5f, 5)]
    public float offsetDist = 1;
    // Start is called before the first frame update
    void Start()
    {
        parent = transform.parent.gameObject;
    }

    // Update is called once per frame
    void Update()
    {
        //position 1 unit from parent in the direction of target
        var offset = target.transform.position - parent.transform.position;
        transform.position = parent.transform.position + offset.normalized * offsetDist;
        transform.LookAt(target.transform.position);
    }
}
