from plotly.subplots import make_subplots
import plotly.graph_objects as go


def plot_result_rates(relaxation_input, fields, result_compilated, drop=None, error=None):
    nfield = len(fields)
    fig = make_subplots(rows=4, cols=nfield, vertical_spacing = 0.05, subplot_titles=("",))
    fig.update_layout(showlegend=True, height=1000, width=1300, hovermode="closest", title_text="Model Free Analysis result", font=dict(size=14), barmode="overlay")
    fig.update_yaxes(exponentformat="power")
    for i, f in enumerate(fields):
        data = relaxation_input[i]
        if drop is not None:
            data = data[~data['num'].isin(drop)]
        data_pred = result_compilated
        if drop is not None:
            data_pred = data_pred[~data_pred['num'].isin(drop)]
        if error is None:
            fig.add_trace(go.Scatter(x=data["num"], y=data["R1"], mode="lines", marker_color="black", name = f"R1 exp {f} MHz", error_y =dict(type="data", array=data["E_R1"], thickness=1), hoverinfo = "x+y+text"), row=1, col=i+1)
            fig.add_trace(go.Scatter(x=data_pred["num"], y=data_pred[f"R1_pred_{f}"],mode="lines", marker_color="red", name = f"R1 pred {f} MHz", hoverinfo = "x+y+text"), row=1, col=i+1)
            fig.add_trace(go.Scatter(x=data["num"], y=data["R2"], mode="lines", marker_color="black", name = f"R2 exp {f} MHz", error_y =dict(type="data", array=data["E_R2"], thickness=1), hoverinfo = "x+y+text"), row=2, col=i+1)
            fig.add_trace(go.Scatter(x=data_pred["num"], y=data_pred[f"R2_pred_{f}"],mode="lines", marker_color="red", name = f"R2 pred {f} MHz", hoverinfo = "x+y+text"), row=2, col=i+1)
            fig.add_trace(go.Bar(x=data_pred["num"], y=data_pred['Rex']*(f**2), marker_color="purple", name = f"Rex pred {f} MHz", hoverinfo = "x+y+text"), row=2, col=i+1)        
            if 'NOE' in data.columns:
                fig.add_trace(go.Scatter(x=data["num"], y=data["NOE"], mode="lines", marker_color="black", name = f"nOe exp {f} MHz", error_y =dict(type="data", array=data["E_NOE"], thickness=1), hoverinfo = "x+y+text"), row=3, col=i+1)
                fig.add_trace(go.Scatter(x=data_pred["num"], y=data_pred[f"NOE_pred_{f}"],mode="lines", marker_color="red", name = f"nOe pred {f} MHz", hoverinfo = "x+y+text"), row=3, col=i+1)
            if 'etaXY' in data.columns:
                fig.add_trace(go.Scatter(x=data["num"], y=data["etaXY"], mode="lines", marker_color="black", name = f"etaXY exp {f} MHz", error_y =dict(type="data", array=data["E_etaXY"], thickness=1), hoverinfo = "x+y+text"), row=4, col=i+1)
                fig.add_trace(go.Scatter(x=data_pred["num"], y=data_pred[f"etaXY_pred_{f}"],mode="lines", marker_color="red", name = f"etaXY pred {f} MHz", hoverinfo = "x+y+text"), row=4, col=i+1)
        if error is not None:
            if drop is not None:
                error = error[~error['num'].isin(drop)]
            fig.add_trace(go.Scatter(x=data["num"], y=data["R1"], mode="lines", marker_color="black", name = f"R1 exp {f} MHz", error_y =dict(type="data", array=data["E_R1"], thickness=1), hoverinfo = "x+y+text"), row=1, col=i+1)
            fig.add_trace(go.Scatter(x=data_pred["num"], y=data_pred[f"R1_pred_{f}"],mode="lines", marker_color="red", name = f"R1 pred {f} MHz", error_y =dict(type="data", array=error[f"R1_pred_{f}"], thickness=1), hoverinfo = "x+y+text"), row=1, col=i+1)
            fig.add_trace(go.Scatter(x=data["num"], y=data["R2"], mode="lines", marker_color="black", name = f"R2 exp {f} MHz", error_y =dict(type="data", array=data["E_R2"], thickness=1), hoverinfo = "x+y+text"), row=2, col=i+1)
            fig.add_trace(go.Scatter(x=data_pred["num"], y=data_pred[f"R2_pred_{f}"],mode="lines", marker_color="red", name = f"R2 pred {f} MHz", error_y =dict(type="data", array=error[f"R2_pred_{f}"], thickness=1), hoverinfo = "x+y+text"), row=2, col=i+1)
            fig.add_trace(go.Bar(x=data_pred["num"], y=data_pred['Rex']*(f**2), marker_color="purple", name = f"Rex pred {f} MHz", hoverinfo = "x+y+text"), row=2, col=i+1)        
            if 'NOE' in data.columns:
                fig.add_trace(go.Scatter(x=data["num"], y=data["NOE"], mode="lines", marker_color="black", name = f"nOe exp {f} MHz", error_y =dict(type="data", array=data["E_NOE"], thickness=1), hoverinfo = "x+y+text"), row=3, col=i+1)
            fig.add_trace(go.Scatter(x=data_pred["num"], y=data_pred[f"NOE_pred_{f}"],mode="lines", marker_color="red", name = f"nOe pred {f} MHz", error_y =dict(type="data", array=error[f"NOE_pred_{f}"], thickness=1), hoverinfo = "x+y+text"), row=3, col=i+1)
            if 'etaXY' in data.columns:
                fig.add_trace(go.Scatter(x=data["num"], y=data["etaXY"], mode="lines", marker_color="black", name = f"etaXY exp {f} MHz", error_y =dict(type="data", array=data["E_etaXY"], thickness=1), hoverinfo = "x+y+text"), row=4, col=i+1)
            fig.add_trace(go.Scatter(x=data_pred["num"], y=data_pred[f"etaXY_pred_{f}"],mode="lines", marker_color="red", name = f"etaXY pred {f} MHz", error_y =dict(type="data", array=error[f"etaXY_pred_{f}"], thickness=1), hoverinfo = "x+y+text"), row=4, col=i+1)     
    fig.update_xaxes(title_text="residue", row=4)
    fig.update_yaxes(title_text="R1 (s\u207b\N{SUPERSCRIPT ONE})", automargin=True, row=1, col=1)
    fig.update_yaxes(title_text="R2 (s\u207b\N{SUPERSCRIPT ONE})", automargin=True, row=2, col=1)
    fig.update_yaxes(title_text="nOe", automargin=False, row=3, col=1)
    fig.update_yaxes(title_text="etaXY (s\u207b\N{SUPERSCRIPT ONE})", automargin=True, row=4, col=1)
    
    return fig

def plot_result_params(result_compilated, drop=None, error=None):
    data = result_compilated
    fig = make_subplots(rows=4, cols=1, vertical_spacing = 0.05, subplot_titles=("",))
    fig.update_layout(showlegend=True, height=2200, width=1300, hovermode="closest", title_text="Model Free Analysis result", font=dict(size=14), barmode="overlay")
    fig.update_yaxes(exponentformat="power")
    if drop is not None:
        data = data[~data['num'].isin(drop)]
    if error is None:
        if 'tm' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["tm"]*1e9, mode="markers+lines", marker_color="blue", name = "tau1", hoverinfo = "x+y+text"), row=1, col=1)
        if 'tf' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["tf"]*1e9, mode="markers+lines", marker_color="red", name = "tau2", hoverinfo = "x+y+text"), row=1, col=1)
        if 'ts' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["ts"]*1e9, mode="markers+lines", marker_color="green", name = "tau3", hoverinfo = "x+y+text"), row=1, col=1)    
        
        if 'S2s' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["S2s"], mode="markers+lines", marker_color="blue", name = "A1", hoverinfo = "x+y+text"), row=2, col=1)
        if 'S2f' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["S2f"], mode="markers+lines", marker_color="red", name = "A2", hoverinfo = "x+y+text"), row=2, col=1)
        if 'S3' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["S3"], mode="markers+lines", marker_color="green", name = "A3", hoverinfo = "x+y+text"), row=2, col=1)

        if 'theta' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["theta"], mode="markers+lines", marker_color="black", name = "theta", hoverinfo = "x+y+text"), row=3, col=1)

        if 'Rex' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["Rex"]*(600**2), mode="markers+lines", marker_color="purple", name = "Rex at 600 Mhz", hoverinfo = "x+y+text"), row=4, col=1)
    else:
        if drop is not None:
            error = error[~error['num'].isin(drop)]
        if 'tm' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["tm"]*1e9, mode="markers+lines", marker_color="red", name = "tau1", error_y =dict(type="data", array=error["E_tm"]*1e9, thickness=1), hoverinfo = "x+y+text"), row=1, col=1)
        if 'tf' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["tf"]*1e9, mode="markers+lines", marker_color="blue", name = "tau2", error_y =dict(type="data", array=error["E_tf"]*1e9, thickness=1), hoverinfo = "x+y+text"), row=1, col=1)
        if 'ts' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["ts"]*1e9, mode="markers+lines", marker_color="green", name = "tau3", error_y =dict(type="data", array=error["E_ts"]*1e9, thickness=1),  visible='legendonly', hoverinfo = "x+y+text"), row=1, col=1)
        
        if 'S2s' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["S2s"], mode="markers+lines", marker_color="red", name = "A1", error_y =dict(type="data", array=error["E_S2s"], thickness=1), hoverinfo = "x+y+text"), row=2, col=1)
        if 'S2f' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["S2f"], mode="markers+lines", marker_color="blue", name = "A2", error_y =dict(type="data", array=error["E_S2f"], thickness=1), hoverinfo = "x+y+text"), row=2, col=1)
        if 'S3' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["S3"], mode="markers+lines", marker_color="green", name = "A3", error_y =dict(type="data", array=error["E_S3"], thickness=1), hoverinfo = "x+y+text"), row=2, col=1)

        if 'theta' in data.columns:
            fig.add_trace(go.Scatter(x=data["num"], y=data["theta"], mode="markers+lines", marker_color="black", name = "theta", error_y =dict(type="data", array=error["E_theta"], thickness=1), hoverinfo = "x+y+text"), row=3, col=1)

        if 'Rex' in data.columns:
            fig.add_trace(go.Bar(x=data["num"], y=data["Rex"]*(600**2), marker_color="purple", name = "Rex at 600 Mhz", error_y =dict(type="data", array=error["E_Rex"]*(600**2), thickness=1), hoverinfo = "x+y+text"), row=4, col=1)        

    fig.update_xaxes(title_text="residue", row=4)
    fig.update_yaxes(title_text="tau (ns\u207b\N{SUPERSCRIPT ONE})", automargin=True, row=1, col=1)
    fig.update_yaxes(title_text="Amplitude", automargin=True, row=2, col=1)
    fig.update_yaxes(title_text="theta (°)", automargin=True, row=3, col=1)
    fig.update_yaxes(title_text="Rex (ns\u207b\N{SUPERSCRIPT ONE})", automargin=True, row=4, col=1)

    return fig