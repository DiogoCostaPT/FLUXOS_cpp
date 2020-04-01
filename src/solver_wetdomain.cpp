
#include <iostream>
#include <armadillo>
#include <memory> 

#include "GlobVar.h"
#include "solver_wetdomain.h"

void solver_wet(GlobVar& ds, unsigned int irow, unsigned int icol){

    unsigned int iw,ie, is,in,inn, NROWSl, NCOLSl,dx,dy;
    double fe1,fe2,fe3,fn1,fn2,fn3,zw,zp,ze,zs,
        zn,znn,hw,he,hs,hn,hnn,fe2p,fn3p,qe,qp,qn,rw,rp,re,rn;
    double dze,dqe,dre,
        hme,qme,rme,dzn;
    double dqn,drn,hmn,hne,qmn,rmn,volrat,ume,vmn,volpot,cf;
    double cme,cmn,vme,umn,//dzbewr,dzbnsr,
        cnp,cne,cnn,up,un,us0,use0,une,vp,ve,vw,vwn,ven,txye,txyn,fe1c,fe2c,fe3c,fn1c,fn2c,fn3c,dzea,dhea;
    double cc,c1,c2,c3,c1a,c2a,c3a,a11,a12,a21,a22,a31,a32,a33,fe1r,fe2r,fe3r,dzna,dhna,a13,a23,fn1r,fn2r,fn3r;
    double hp ,hp0; 
    double zbw,zbp,zbe,zbs,zbn,zbnn,zbpe,zbpn; //zbpw,zbps
    double dtl, hdryl,gaccl,kspl, nueml, cvdefl; 
    float ldp,lde,ldn;
    bool lroe;

    NROWSl = ds.NROWS;
    NCOLSl = ds.NCOLS;
    hdryl = ds.hdry;
    gaccl = ds.gacc;
    kspl = (*ds.ks).at(irow,icol);
    nueml = ds.nuem;
    cvdefl = ds.cvdef;
    
    is=icol-1;
    in=icol+1;
    inn=fmin(icol+2,NCOLSl+1);
    lroe = true;
    iw=irow-1;
    ie=irow+1;
   
    dx  = ds.dxy;
    dy  = ds.dxy;
    
    ldp = (*ds.ldry).at(irow,icol);
    lde = (*ds.ldry).at(ie,icol);
    ldn = (*ds.ldry).at(irow,in);

    // CELL CENTER VALUES
    zbw = (*ds.zb).at(iw,icol);
    zbp = (*ds.zb).at(irow,icol);
    zbe = (*ds.zb).at(ie,icol);
    zbs = (*ds.zb).at(irow,is);
    zbn = (*ds.zb).at(irow,in);
    zbnn= (*ds.zb).at(irow,inn);
    zw=(*ds.z).at(iw,icol);
    zp=(*ds.z).at(irow,icol);
    ze=(*ds.z).at(ie,icol);
    zs=(*ds.z).at(irow,is);
    zn=(*ds.z).at(irow,in);
    znn=(*ds.z).at(irow,inn);
    qp=(*ds.qx).at(irow,icol);
    qe=(*ds.qx).at(ie,icol);
    qn=(*ds.qx).at(irow,in);
    rw=(*ds.qy).at(iw,icol);
    rp=(*ds.qy).at(irow,icol);
    re=(*ds.qy).at(ie,icol);
    rn=(*ds.qy).at(irow,in);

    // zbpw=.5*(zbw+zbp);
    zbpe=.5*(zbe+zbp);
    // zbps=.5*(zbs+zbp);
    zbpn=.5*(zbn+zbp);
    hp  = std::max(0.0,(*ds.z).at(irow,icol)-zbp);
    hp0 = std::max(std::max(hdryl,hp),kspl);
    hw=std::max(0.0,zw-zbw);
    he=std::max(0.0,ze-zbe);
    hs=std::max(0.0,zs-zbs);
    hn=std::max(0.0,zn-zbn);
    hnn=std::max(0.0,znn-zbnn);

    // CELL FACE VALUES                    
    hme=.5*(hp+he);          
    qme=.5*(qp+qe);
    rme=.5*(rp+re);
    hmn=.5*(hp+hn);   
    qmn=.5*(qp+qn);
    rmn=.5*(rp+rn);

    dze=ze-zp;
    dqe=qe-qp;
    dre=re-rp;
    dzn=zn-zp;
    drn=rn-rp;
    dqn=qn-qp;
    
    // TIMESTEP
    dtl=ds.dtfl;

    // CELLS WITH SOME DRY NEIGHBOURS     
    if(lde==1.0f) 
    { 
        hme=std::max(0.0,zp-zbpe);
        he=0.0f;;
        qe=0.0f;
        re=0.0f;
        if(hme>hdryl) 
        {
             if(ze>=zp) {
                qme=0.0f;
                rme=0.0f;
            }else 
            {
                dze=fmin(fabs(dze),hme);
                qme=0.296*dze*sqrt(gaccl*dze);    // Ritter solution
                volrat=0.5*dx*hp/dtl;            // available volume rate per m of cell P [m2/s]
                qme=fmin(qme,volrat);
                rme=0.0f;
            }
        }else 
        {
            qme=0.0f;
            rme=0.0f;
            hme=0.0f;            
        }
        lroe= false;
    }
    if(ldn==1.0f) 
    {
        hmn=std::max(0.0,zp-zbpn);
        hn=0.0f;
        qn=0.0f;
        rn=0.0f;
        if(hmn>hdryl) 
        {
            if(zn>=zp) 
            {
                qmn=0.0f;
                rmn=0.0f;
            }else 
            {
                dzn=fmin(fabs(dzn),hmn);
                rmn=0.296*dzn*sqrt(gaccl*dzn);    // Ritter solution
                qmn=0.0f;
                volrat=0.5*dy*hp/dtl;            // available volume rate per m of cell P [m2/s]
                rmn=fmin(rmn,volrat);
            }
        }else 
        {
            qmn=0.0f;
            rmn=0.0f;
            hmn=0.0f;
        }
        lroe=false;
    }

    // CALC TURBULENT STRESS
    cme=sqrt(gaccl*hme);
    cmn=sqrt(gaccl*hmn);
    ume=qme/std::max(hme,hdryl);
    vme=rme/std::max(hme,hdryl);
    umn=qmn/std::max(hmn,hdryl);
    vmn=rmn/std::max(hmn,hdryl);
    
    cnp=cvdefl*(*ds.us).at(irow,icol)*hp+nueml;
    cne=cvdefl*(*ds.us).at(ie,icol)*he+nueml;
    cnn=cvdefl*(*ds.us).at(irow,in)*hn+nueml;
    hne=.5*(cnp+cne)*sqrt(hp*he);
    hnn=.5*(cnp+cnn)*sqrt(hp*hn);

    up=qp/hp0;
    un=qn/std::max(std::max(hn,hdryl),(*ds.ks).at(irow,in));
    us0=(*ds.qx).at(irow,is)/std::max(std::max(hs,hdryl),(*ds.ks).at(irow,is));
    use0=(*ds.qx).at(ie,is)/std::max(std::max((*ds.h).at(ie,is),hdryl),(*ds.ks).at(ie,is));
    une=(*ds.qx).at(ie,in)/std::max(std::max((*ds.h).at(ie,in),hdryl),(*ds.ks).at(ie,in));
    vp=rp/hp0;
    ve=re/std::max(std::max(he,hdryl),(*ds.ks).at(ie,icol));
    vw=rw/std::max(std::max(hw,hdryl),(*ds.ks).at(iw,icol));
    vwn=(*ds.qy).at(iw,in)/std::max(std::max((*ds.h).at(iw,in),hdryl),(*ds.ks).at(iw,in));
    ven=(*ds.qy).at(ie,in)/std::max(std::max((*ds.h).at(ie,in),hdryl),(*ds.ks).at(ie,in));

    txye=hne*((ve-vp)/fabs(dx)+.25*(un+une-us0-use0)/dy);
    txyn=hnn*((un-up)/fabs(dy)+.25*(ve+ven-vw-vwn)/dx);

    // CALC OF CONVECTION FLUXES
    fe1c=qme;
    fe2c=qme*ume;
    fe3c=qme*vme -txye;
    fn1c=rmn;
    fn2c=rmn*umn -txyn;
    fn3c=rmn*vmn;

    // ROE's DISSIPASSION
    if(lroe) 
    { 
        if(hme>ds.hdry) 
        {
            dzea=fabs(dze);
            if(dzea>0.5*hme) 
            {
                dhea=fabs(he-hp);
                dzea=fmin(dzea,dhea);
                dze=std::copysign(dzea,dze);
            }
            cc=.25/cme;
        }else 
        {
            cc=0.0f;
        }
        
        c1=ume;
        c2=ume+cme;
        c3=ume-cme;
        c1a=fabs(c1);
        c2a=fabs(c2);
        c3a=fabs(c3);
        a11=c2*c3a-c2a*c3;
        a12=c2a-c3a;
        a21=c2*c3*(c3a-c2a);
        a22=c2a*c2-c3a*c3;
        a31=vme*(c2*c3a-2.*cme*c1a-c2a*c3);
        a32=vme*(c2a-c3a);
        a33=2.*cme*c1a;

        fe1r=-(a11*dze+a12*dqe)*cc;
        fe2r=-(a21*dze+a22*dqe)*cc;
        fe3r=-(a31*dze+a32*dqe+a33*dre)*cc;

        if(ldp==0.0f&&hmn>ds.hdry) 
        {
                dzna=fabs(dzn);
                dhna=fabs(hn-hp);
                dzna=fmin(dzna,dhna);
                dzn=std::copysign(dzna,dzn);
                cc=.25/cmn;
        }else 
        {
            cc=0.0f;
        };
        c1=vmn;
        c2=vmn+cmn;
        c3=vmn-cmn;
        c1a=std::fabs(c1);
        c2a=std::fabs(c2);
        c3a=std::fabs(c3);
        a11=c2*c3a-c2a*c3;
        a13=c2a-c3a;
        a21=umn*(c2*c3a-2.*cmn*c1a-c2a*c3);
        a22=2.*cmn*c1a;
        a23=umn*(c2a-c3a);
        a31=c2*c3*(c3a-c2a);
        a33=c2a*c2-c3a*c3;

        fn1r=-(a11*dzn+a13*drn)*cc;
        fn2r=-(a21*dzn+a22*dqn+a23*drn)*cc;
        fn3r=-(a31*dzn+a33*drn)*cc;
    } else
    {
        fe1r=0.0f;
        fe2r=0.0f;
        fe3r=0.0f;
        fn1r=0.0f;
        fn2r=0.0f;
        fn3r=0.0f;
    }
    
    // PRESSURE AT CELL SIDE
    fe2p=.5*gaccl*(hme*hme);
    fn3p=.5*gaccl*(hmn*hmn); 

    // SUM OF ALL FLUXES
    fe1=fe1c+fe1r;
    fe2=fe2c+fe2r+fe2p;
    fe3=fe3c+fe3r;
    fn1=fn1c+fn1r;
    fn2=fn2c+fn2r;
    fn3=fn3c+fn3r+fn3p;
            
    // BOUNDARY CONDITIONS (WEIR DISCHARGE RATE) 
    if (icol==1 || icol==NCOLSl)
    {
        fn1=std::min(volrat,sqrt(gaccl)*pow(std::fmax(hp,0.0f),1.5));
    }
    if (irow==1 || irow==NROWSl)
    {
        fe1=std::min(volrat,sqrt(gaccl)*pow(std::fmax(hp,0.0f),1.5));
    }

    // CHECK MASS BALANCE (restrict outflow flux to available water)        
    volpot=ds.arbase*hp;    // volume in cell P [m3]
    volrat=volpot/dtl;      // max flux rate 
        
    if(volrat>0.0f) { // cell has water
        if(fe1>0.0f&&fn1>0.0f) 
        {
            if(fe1*dy+fn1*dx>volrat) {
                cf=fn1*dx/(fe1*dy+fn1*dx);
                fe1=(1.-cf)*volrat/dy;            // [m3/s]
                fn1=cf*volrat/dx;
            }
        }else if(fe1>0.0f) 
        {
            fe1=fmin(fe1*dy,volrat)/dy; 
        }else if(fn1>0.0f) 
        {
            fn1=fmin(fn1*dx,volrat)/dx;
        }
    }else // cell has no water 
    {
        fe1=0.0f;
        fn1=0.0f;
    }

    // SAVE MASS AND MOMENTUM FLUXES
    (*ds.fn_1).at(irow,icol)=fn1;
    (*ds.fn_2).at(irow,icol)=fn2;
    (*ds.fn_3).at(irow,icol)=fn3;
    (*ds.fe_1).at(irow,icol)=fe1;
    (*ds.fe_2).at(irow,icol)=fe2;
    (*ds.fe_3).at(irow,icol)=fe3;
    
} 

