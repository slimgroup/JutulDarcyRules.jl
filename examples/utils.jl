#### Patchy saturation model

function Patchy(sw::AbstractMatrix{T}, vp::AbstractMatrix{T}, rho::AbstractMatrix{T}, phi::AbstractMatrix{T};
    bulk_min::T=T(36.6f9), bulk_fl1::T=T(2.735f9), bulk_fl2::T=T(0.125f9), ρw::T=T(501.9f0), ρo::T=T(1053.0f0)) where T

    ### works for channel problem
    vs = vp./sqrt(3f0)
    bulk_sat1 = rho .* (vp.^2f0 - 4f0/3f0 .* vs.^2f0)
    shear_sat1 = rho .* (vs.^2f0)

    patch_temp = bulk_sat1 ./(bulk_min .- bulk_sat1) - 
    bulk_fl1 ./ phi ./ (bulk_min .- bulk_fl1) + 
    bulk_fl2 ./ phi ./ (bulk_min .- bulk_fl2)

    bulk_sat2 = bulk_min./(1f0./patch_temp .+ 1f0)

    bulk_new = 1f0./( (1f0.-sw)./(bulk_sat1+4f0/3f0*shear_sat1) 
    + sw./(bulk_sat2+4f0/3f0*shear_sat1) ) - 4f0/3f0*shear_sat1

    rho_new = rho + phi .* sw * (ρw - ρo)

    Vp_new = sqrt.((bulk_new+4f0/3f0*shear_sat1)./rho_new)
    return Vp_new, rho_new

end

function Patchy(sw::AbstractArray{T, 3}, vp::AbstractMatrix{T}, rho::AbstractMatrix{T}, phi::AbstractMatrix{T};
    bulk_min::T=T(36.6f9), bulk_fl1::T=T(2.735f9), bulk_fl2::T=T(0.125f9), ρw::T=T(501.9f0), ρo::T=T(1053.0f0)) where T

    stack = [Patchy(sw[i,:,:], vp, rho, phi; bulk_min=bulk_min, bulk_fl1=bulk_fl1, bulk_fl2=bulk_fl2, ρw = ρw, ρo=ρo) for i = 1:size(sw,1)]
    return [stack[i][1] for i = 1:size(sw,1)], [stack[i][2] for i = 1:size(sw,1)]
end