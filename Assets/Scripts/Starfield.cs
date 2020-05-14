using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Assertions;


public class Starfield : MonoBehaviour
{
    public int maxStars = 100;
    public float starFieldSize = 10;

    public float starSize = 1f;
    public float starSizeRange = 5f;
    public bool colorize = false;

    ParticleSystem particles;
    ParticleSystem.Particle[] stars;
    private ParticleSystem.EmitParams emitParams;

    void Awake()
    {
        stars = new ParticleSystem.Particle[maxStars];
        particles = GetComponent<ParticleSystem>();

        Assert.IsNotNull(particles, "Particle system missing from object!");

        for (int i = 0; i < maxStars; i++)
        {
            float randSize = Random.Range(1, starSizeRange );                       // Randomize star size within parameters
            float scaledColor = (true == colorize) ? randSize - starSizeRange : 1f;         // If coloration is desired, color based on size

            stars[i].position = Random.insideUnitSphere.normalized * starFieldSize;
            stars[i].startSize = starSize * randSize;
            stars[i].startColor = new Color(1f, scaledColor, scaledColor, 1f);
        }
        particles.SetParticles(stars, stars.Length);

        emitParams.position = new Vector3(0, 2, 0);
        emitParams.velocity = new Vector3(0, 0, 0);
        emitParams.startLifetime = float.MaxValue;

        //particles.Emit(emitParams, stars.Length);
        //particles.Play();
    }

}
