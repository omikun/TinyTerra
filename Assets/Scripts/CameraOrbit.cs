// from https://www.youtube.com/watch?v=bVo0YLLO43s
// Emergent Sagas
// Simple Mouse Orbit Camera Controller How To - Quick Unity Tutorials
// https://pastebin.com/raw/2RX8fpJ3
using UnityEngine;
using System.Collections;

public class CameraOrbit : MonoBehaviour 
{

    protected Transform _XForm_Camera;
    protected Transform parent;
    private Vector3 startOffset;
    private Quaternion startParentRotation;
    private Quaternion startRotation;
    [Tooltip("Speed at which the camera rotates. (Camera uses Slerp for rotation.)")]
    public float rotateSpeed = 90.0f;

    protected Vector3 _LocalRotation;
    protected float _CameraDistance = 10f;

    public float MouseSensitivity = 4f;
    public float ScrollSensitvity = 2f;
    public float OrbitDampening = 10f;
    public float ScrollDampening = 6f;

    public bool OrbitCameraEnabled = false;


    // Use this for initialization
    void Start() {
        this._XForm_Camera = this.transform;
        this.parent = this.transform.parent;

        //start in lag camera mode
        startOffset = transform.localPosition;
        startRotation = transform.localRotation;
        startParentRotation = parent.localRotation;
        transform.SetParent(null);
    }


    void LateUpdate() {
        //OrbitCameraEnabled = Input.GetMouseButton(0);
        bool prevMode = OrbitCameraEnabled;
        OrbitCameraEnabled ^= Input.GetKeyDown(KeyCode.Space);

        if (!prevMode && OrbitCameraEnabled)
        {
            transform.SetParent(parent);
        }
        if (prevMode && !OrbitCameraEnabled)
        {
            transform.SetParent(null);
            transform.localPosition = startOffset;
            transform.localRotation = startRotation;
            parent.localRotation = startParentRotation;
        }

        if (OrbitCameraEnabled)
        { 
            OrbitCamera();
        }
        else
        {
            LagCamera();
        }
    }
    
    void LagCamera() 
    {
        if (parent != null)
        {
            transform.position = parent.TransformPoint(startOffset);
            transform.rotation = Quaternion.Slerp(transform.rotation, parent.rotation, rotateSpeed * Time.deltaTime);
        }
        
    }

    void OrbitCamera()
    {
        //Rotation of the Camera based on Mouse Coordinates
        if (Input.GetMouseButton(0) && (Input.GetAxis("Mouse X") != 0 || Input.GetAxis("Mouse Y") != 0))
            {
                _LocalRotation.x += Input.GetAxis("Mouse X") * MouseSensitivity;
                _LocalRotation.y += -1 * Input.GetAxis("Mouse Y") * MouseSensitivity;

                //Clamp the y Rotation to horizon and not flipping over at the top
                if (_LocalRotation.y < -89f)
                    _LocalRotation.y = -89f;
                else if (_LocalRotation.y > 90f)
                    _LocalRotation.y = 90f;
            }
        //Zooming Input from our Mouse Scroll Wheel
        if (Input.GetAxis("Mouse ScrollWheel") != 0f)
        {
            float ScrollAmount = Input.GetAxis("Mouse ScrollWheel") * ScrollSensitvity;

            ScrollAmount *= (this._CameraDistance * 0.3f);

            this._CameraDistance += ScrollAmount * -1f;

            this._CameraDistance = Mathf.Clamp(this._CameraDistance, 1.5f, 100f);
        }

        //Actual Camera Rig Transformations
        Quaternion QT = Quaternion.Euler(_LocalRotation.y, _LocalRotation.x, 0);
        this.parent.rotation = Quaternion.Lerp(this.parent.rotation, QT, Time.deltaTime * OrbitDampening);

        if ( this._XForm_Camera.localPosition.z != this._CameraDistance * -1f )
        {
            this._XForm_Camera.localPosition = new Vector3(0f, 0f, Mathf.Lerp(this._XForm_Camera.localPosition.z, this._CameraDistance * -1f, Time.deltaTime * ScrollDampening));
        }

    }
}