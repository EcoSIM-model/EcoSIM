"""Check current-hour temperature and single advancement in hourly crop history.

Use an emerged, pre-flowering, single-branch window with variable temperatures.
The history origin must be supplied explicitly because these files need not
contain a CF time coordinate. No model or parameter files are modified.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from netCDF4 import Dataset


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history',type=Path,required=True)
    parser.add_argument('--origin',required=True)
    parser.add_argument('--start',required=True)
    parser.add_argument('--end',required=True,help='Exclusive interval end')
    parser.add_argument('--offset',type=float,required=True)
    parser.add_argument('--xrni',type=float,required=True)
    parser.add_argument('--xrla',type=float,required=True)
    parser.add_argument('--pft',type=int,default=0)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    with Dataset(args.history) as ds:
        bounds=np.asarray(ds['time_bounds'][:],dtype=float)
        assert np.allclose(np.diff(bounds,axis=1),1/24,atol=1e-9)
        assert np.allclose(bounds[1:,0],bounds[:-1,1],atol=1e-9)
        ix=pd.Timestamp(args.origin)+pd.to_timedelta(np.rint(bounds[:,0]*24).astype(int),unit='h')
        names=['TC_Groth_pft','CanPhenol_TSTRESS_pft','MainBranchNodeNumber_pft',
               'ShootNodeNumber_pft','Num_Leaves_pft','RNodeAppear_pft','RLeafAppear_pft',
               'OXY_STRESS_pft','PSI_CAN_pft','Growth_Stage_pft']
        h=pd.DataFrame({k:np.ma.filled(ds[k][:,args.pft],np.nan).astype(float) for k in names},index=ix)
    keep=(h.index>=pd.Timestamp(args.start))&(h.index<pd.Timestamp(args.end))
    z=h.loc[keep]
    assert len(z)>24 and np.isfinite(z.to_numpy()).all()
    assert z.Growth_Stage_pft.eq(1).all(),'Require emerged, pre-flowering crop'
    assert np.max(abs(z.MainBranchNodeNumber_pft-z.ShootNodeNumber_pft))<1e-6,'Require single-branch window'
    def response(tc):
        tk=tc+273.15+args.offset;rt=8.3143*tk;st=710*tk
        return np.exp(24.269-60000/rt)/(1+np.exp((197500-st)/rt)+np.exp((st-218500)/rt))
    error=z.CanPhenol_TSTRESS_pft-response(z.TC_Groth_pft)
    old_error=z.CanPhenol_TSTRESS_pft-response(h.TC_Groth_pft.shift().loc[keep])
    # Growth temperature is truncated to .001 C and multipliers are float32.
    assert error.abs().max()<1e-4,('current temperature mismatch',error.abs().max())
    assert old_error.abs().mean()>1e-3,'Window must distinguish current from retained temperature'
    modifier=np.exp(.025*np.maximum(z.PSI_CAN_pft,-1000))*np.sqrt(z.OXY_STRESS_pft)
    expected_node=args.xrni*z.CanPhenol_TSTRESS_pft*modifier
    expected_leaf=args.xrla*z.CanPhenol_TSTRESS_pft*modifier
    node_delta=h.MainBranchNodeNumber_pft.diff().loc[keep]
    leaf_delta=h.Num_Leaves_pft.diff().loc[keep]
    # Absolute counter output precision limits hourly differences, especially
    # once several nodes have accumulated. A duplicated update is much larger.
    node_error=(node_delta-expected_node).abs().max()
    leaf_error=(leaf_delta-expected_leaf).abs().max()
    assert node_error<2e-6,('node advancement mismatch',node_error)
    assert leaf_error<2e-6,('leaf advancement mismatch',leaf_error)
    assert np.max(abs(z.RNodeAppear_pft/24-expected_node))<1e-7
    assert np.max(abs(z.RLeafAppear_pft/24-expected_leaf))<1e-7
    report=dict(status='passed',hours=len(z),origin=args.origin,start=args.start,end_exclusive=args.end,
                offset_C=args.offset,current_temperature_max_abs_error=float(error.abs().max()),
                previous_temperature_mean_abs_error=float(old_error.abs().mean()),
                node_increment_max_abs_error=float(node_error),leaf_increment_max_abs_error=float(leaf_error),
                same_hour_temperature=True,single_hourly_node_and_leaf_increment=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
