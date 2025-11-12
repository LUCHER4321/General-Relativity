using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GeodesicMotion : MonoBehaviour
{
    public GeneralRelativity relativityConfig;
    public FourVector position;
    public FourVector velocity;

    //d²x^λ/dτ² + Γ^λ_μν × (dx^μ/dτ) × (dx^ν/dτ) = 0
    private void IntegrateGeodesic()
    {
        SpaceTimeMetric metric = MetricCalculator.CalculateMetric(position, relativityConfig);
        FourVector acceleration = new FourVector(0, 0, 0, 0);
        for (int lambda = 0; lambda < 4; lambda++)
        {
            double sum = 0;
            for (int mu = 0; mu < 4; mu++)
            {
                for (int nu = 0; nu < 4; nu++)
                {
                    sum -= metric.christoffelSymbols[lambda, mu, nu] *
                           velocity.components[mu] * velocity.components[nu];
                }
            }
            acceleration.components[lambda] = sum;
        }
        for (int i = 0; i < 4; i++)
        {
            velocity.components[i] += acceleration.components[i] * Time.fixedDeltaTime;
            position.components[i] += velocity.components[i] * Time.fixedDeltaTime;
        }
    }

    void FixedUpdate()
    {

    }
}
