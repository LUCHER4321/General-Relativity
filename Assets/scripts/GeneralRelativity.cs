using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public enum MetricType
{
    Minkowski,
    Schwarzschild,
    Kerr
}

public class GeneralRelativity : MonoBehaviour
{
    public static double lightSpeed = 299792458.0; // c (m/s)
    public static double c2 = Math.Pow(lightSpeed, 2.0);
    public static double gravitationalConstant = 6.67430e-11; // G (m³/(kg×s²))

    public MetricType currentMetric = MetricType.Minkowski;
    public Rigidbody rb;

    public double GetKerrA()
    {
        Matrix4x4 rotationMatrix = Matrix4x4.Rotate(rb.inertiaTensorRotation);
        Matrix4x4 scaleMatrix = Matrix4x4.Scale(rb.inertiaTensor);
        Matrix4x4 inertiaTensor = rotationMatrix * scaleMatrix * rotationMatrix.inverse;
        Vector3 angularMomentum = inertiaTensor.MultiplyVector(rb.angularVelocity);
        return angularMomentum.magnitude / lightSpeed;
    }
}
