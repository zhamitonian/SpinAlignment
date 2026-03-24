#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TVirtualFitter.h"

using namespace std;

void cos_theta() {
    gStyle->SetOptStat(0); // 关闭统计信息

    double efficiency[13][10] = {0};
    double efficiencyError[13][10] = {0};
    double nsig[13][10] = {0};
    double nsigError[13][10] = {0};
    double nsig_real[13][10] = {0};
    double nsigError_real[13][10] = {0};
    
    
    // 打开输入文件
    ifstream file1("/besfs5/groups/tauqcd/luruihua/cross_check/draw/efficiency/nsig_truth_wangz.txt");
    //ifstream file1("/besfs5/groups/tauqcd/luruihua/cross_check/fit/data_try/nsig_results.txt");
    ifstream file2("/besfs5/groups/tauqcd/luruihua/cross_check/draw/efficiency/efficiency_results.txt");

    if (!file1.is_open()) {
        cerr << "Error opening nsig_results!" << endl;
        return;
    }
    if (!file2.is_open()) {
        cerr << "Error opening efficiency_results!" << endl;
        return;
    }

    // 读取前两个文件生成第一个数据集
    string line;
    int row1 = 0;
    int x, y;
    double max[13]= {0},max1[13]= {0};
    while (getline(file2, line) && row1 < 130) {
        stringstream ss2(line);
        double c, d, e, f;
        ss2 >> c >> d >> e >> f;

        getline(file1, line);
        stringstream ss1(line);
        double a, b, g, h;
        ss1 >> a >> b >> g >> h;

        x = floor(row1 / 10);
        y = row1 - floor(row1 / 10) * 10;

        efficiency[x][y] = e ;
        if(e == 0) efficiency[x][y] = 1; // 避免除以0
        efficiencyError[x][y] = f ;

        nsig[x][y] = round(g);
        if(nsig[x][y]<0) {nsig[x][y]=0;nsigError[x][y]=0;}
        nsigError[x][y] = round(h*10.0)/10.0;

        /*
        nsig_real[x][y] = nsig[x][y] / efficiency[x][y];
        nsigError_real[x][y] = sqrt(nsigError[x][y] * nsigError[x][y] / (efficiency[x][y] * efficiency[x][y]) +
                                    nsig[x][y] * nsig[x][y] * efficiencyError[x][y] * efficiencyError[x][y] /
                                    (efficiency[x][y] * efficiency[x][y] * efficiency[x][y] * efficiency[x][y]));
        */
        nsig_real[x][y] = nsig[x][y] ;
        nsigError_real[x][y] = nsigError[x][y] ;
        nsig_real[x][y] = round(nsig_real[x][y]);
        nsigError_real[x][y] = round(nsigError_real[x][y]*10.0)/10.0;

        if(nsig[x][y] > max1[x]) {
            max1[x] = nsig[x][y];
        }
        if(nsig_real[x][y] > max[x]) {
            max[x] = nsig_real[x][y];
        }
        //if(x==0) {max[x] = 10000;} // 特例处理，避免除以0
        row1++;
    }

    std::ofstream outfile("nsig_real_results.txt");
    if (outfile.is_open()) {
        for (int i = 0; i < 13; ++i) {
            for (int j = 0; j < 10; ++j) {
                outfile << (i + 1) << " " << j <<
                " " << nsig_real[i][j] << " " << nsigError_real[i][j] << std::endl;
            }
        }
        outfile.close();

    } else {
        std::cerr << "Error: Cannot open output file nsig_real_results.txt" << std::endl;
    }

    double x1[10],x_err[10]={0};
        for (int j = 0; j < 10; j++) {
            x1[j] = -0.9+0.2*j; // cos{theta^{star}} 的计算
            x_err[j]=0.0;
        }
    
        double momentum[13] = {0};
        double momentum_center[13] = {0};
        double rho00_values[13] = {0};
        double rho00_errors[13] = {0};

    double num = 0, num_error;
    TGraphErrors* h[13];
    for (int i = 0; i < 13; i++) {
        /*
        if(i<=2) h[i] = new TGraphErrors(4, x1+3, nsig_real[i]+3, x_err, nsigError_real[i]+3 );// 创建图形对象     
        else
        if(i>=2&&i<=5) h[i] = new TGraphErrors(6, x1+2, nsig_real[i]+2, x_err, nsigError_real[i]+2 );// 创建图形对象
        else
        if(i>5) h[i] = new TGraphErrors(8, x1+1, nsig_real[i]+1, x_err, nsigError_real[i]+1 );// 创建图形对象     
        else
            
        if(i==1) 
        h[i] = new TGraphErrors(6, x1+2, nsig_real[i]+2, x_err, nsigError_real[i]+2 );// 创建图形对象
        else
        
        if(i<=7) 
        h[i] = new TGraphErrors(8, x1+1, nsig_real[i]+1, x_err+1, nsigError_real[i]+1 );// 创建图形对象     
        else
        */
        h[i] = new TGraphErrors(10, x1, nsig_real[i], x_err, nsigError_real[i] );// 创建图形对象     
        TCanvas* c = new TCanvas("c", "c", 800, 600);
        c->cd();

        for (int j = 0; j < 10; j++) {
            cout << "Tree Group: " << (i + 1) << ", Tree Index: " << j 
                 <<" Nsig:" << nsig[i][j]
                 << ", eff: " << efficiency[i][j]
                 << ", Nsig Real: " << nsig_real[i][j] 
                 << ", Nsig Real Error: " << nsigError_real[i][j]
                 << endl;
        }
        
        
        // 设置坐标轴标题和标度
        h[i]->SetTitle(Form("%.1f<P(Ks)<%.1f", (i + 1)/10.0, (i + 2)/10.0));
        h[i]->GetXaxis()->SetTitle("cos#theta^{*}");
        h[i]->GetYaxis()->SetTitle("N_{obs}/#epsilon");
        h[i]->GetXaxis()->SetTitleOffset(1.2); // 调整X轴标题偏移
        h[i]->GetYaxis()->SetTitleOffset(1.2); // 调整Y轴标题偏移
        h[i]->GetXaxis()->CenterTitle(); // X轴标题居中
        h[i]->GetYaxis()->CenterTitle(); // Y轴标题居中
        h[i]->GetXaxis()->SetTitleSize(0.04); // 设置X轴标题大小
        h[i]->GetYaxis()->SetTitleSize(0.04); // 设置Y轴标题大小
        h[i]->GetXaxis()->SetLabelSize(0.03); // 设置X轴标签大小
        h[i]->GetYaxis()->SetLabelSize(0.03); // 设置Y轴标签大小
        h[i]->SetMinimum(0); // 设置Y轴最小值为0
        /*
        if(i==0) h[i]->SetMaximum(12000); 
        else
        if(i==1) h[i]->SetMaximum(30000); 
        else
        */
        h[i]->SetMaximum(1.8 * max[i]); 
        h[i]->GetXaxis()->SetLimits(-1, 1); // 设置X轴范围

        // 设置样式
        h[i]->SetLineColor(kBlack);
        h[i]->SetMarkerStyle(20);
        h[i]->SetMarkerColor(kBlack);
        h[i]->SetMarkerSize(0.8); // 设置数据点大小
        h[i]->Draw("AP"); // E1显示误差条，SAME表示叠加绘制

        TF1 func("fitFunc", "[1]*0.75*((1-[0]+(3*[0]-1)*x*x))", -1, 1);
        func.SetParameter(0, 0.5); // 设置初始参数值
        func.SetParameter(1, 15000); // 设置初始参数值
        TFitResultPtr fitResult = h[i]->Fit(&func, "S", "", -1, 1);
        double rho00 = fitResult->Parameter(0);
        double rho00Error = fitResult->ParError(0);
        double rho00Chi2 = fitResult->Chi2();

        momentum[i] = i + 1;  // 动量bin编号
        momentum_center[i] = (i + 1.5) / 10.0;  // 动量bin中心值 (GeV/c)
        rho00_values[i] = rho00;
        rho00_errors[i] = rho00Error;

        // 设置图例
        TLegend* legend = new TLegend(0.55, 0.65, 0.85, 0.85,"","BRNDC"); // 图例位置
        double nBins = h[i]->GetN();
        double nParams = fitResult->NPar();
        double ndf = nBins - nParams;
        legend->AddEntry((TObject*)0, Form("#rho_{00} = %.4f #pm %.4f", rho00, rho00Error), "");
        legend->AddEntry((TObject*)0, Form("#chi^{2}/ndf = %.1f / %.0f = %.1f", rho00Chi2, ndf ,rho00Chi2 / ndf), "");
        legend->SetFillColor(10); // 设置图例背景
        legend->SetTextSize(0.035); // 设置图例文本大小
        legend->SetBorderSize(0); // 去掉图例边框
        //legend->Draw("same"); 
        
        int n_points = 100;
        TGraphErrors* error_band = new TGraphErrors(n_points);
        double x_min = h[i]->GetXaxis()->GetXmin();
        double x_max = h[i]->GetXaxis()->GetXmax();
        cout << "x_min: " << x_min << ", x_max: " << x_max << endl;
        for (int j = 0; j < n_points; j++) {
            double x = x_min + j * (x_max - x_min) / (n_points - 1);
            //cout << "Point " << j << ": x = " << x << endl;
            error_band->SetPoint(j ,x ,0);
        }

        TVirtualFitter* fitter = TVirtualFitter::GetFitter();
        if (fitter) {
            fitter->GetConfidenceIntervals(error_band,0.683); // 1 sigma
        }
        error_band->SetFillColorAlpha(kBlue - 7, 0.35);
        error_band->SetLineColor(kBlue - 7);
        //error_band->Draw("3 same");
        //legend->AddEntry(error_band, "1 #sigma error band", "f");
        legend->Draw();

        c->SaveAs(Form("cos_theta_%d.pdf", i +1));
        delete legend; // 清理图例资源
        //delete c; // 清理画布资源
    
    }

    TH2F* hist = new TH2F("nsig_real", "real N_{sig}(Ks);P(Ks)GeV/c;cos#theta*", 13, 0.1, 1.4, 10, -1, 1);
    for (int i = 0; i < 13; ++i) {
        for (int j = 0; j < 10; ++j) {
            hist->SetBinContent(i + 1, j + 1, nsig_real[i][j]);
            hist->SetBinError(i + 1, j + 1, nsigError_real[i][j]);
        }
    }
    TCanvas* c1 = new TCanvas("c1", "c1", 1200, 900);
    hist->SetContour(100); 
    hist->GetYaxis()->SetTitleOffset(1.2); 
    hist->Draw("COLZ");
    for (int i = 0; i < 13; ++i) {
        for (int j = 0; j < 10; ++j) {
            double xCenter = hist->GetXaxis()->GetBinCenter(i + 1);
            double yCenter = hist->GetYaxis()->GetBinCenter(j + 1);
            char text[50];
            if(nsig_real[i][j]==0) continue;
            snprintf(text, sizeof(text), "%.0f", nsig_real[i][j]);
            TText* t = new TText(xCenter, yCenter, text);
            t->SetTextSize(0.025); // 设置文本大小
            t->SetTextAlign(22); // 文本居中对齐
            t->Draw();
        }
    }
    c1->SaveAs("nsig_real_data.png");

    TH2F *hist_data = new TH2F("nsig_data", "Data N_{sig}(Ks);P(Ks)GeV/c;cos#theta*", 13, 0.1, 1.4, 10, -1, 1);
    for (int i = 0; i < 13; ++i) {
        for (int j = 0; j < 10; ++j) {          
            if(nsig[i][j]==0) continue;
            hist_data->SetBinContent(i + 1, j + 1, nsig[i][j]);
        }
    }
    TCanvas* canvas_data = new TCanvas("canvas_data", "Data Nsig 2D Distribution", 1200, 900);
    hist_data->SetContour(100);
    hist_data->GetYaxis()->SetTitleOffset(1.2);
    hist_data->Draw("COLZ");
    for (int i = 0; i < 13; ++i) {
        for (int j = 0; j < 10; ++j) {
            double xCenter = hist_data->GetXaxis()->GetBinCenter(i + 1);
            double yCenter = hist_data->GetYaxis()->GetBinCenter(j + 1);
            char text[50];
            if(nsig[i][j]==0) continue;
            snprintf(text, sizeof(text), "%.0f", nsig[i][j]);
            TText* t = new TText(xCenter, yCenter, text);
            t->SetTextSize(0.025); // 设置文本大小
            t->SetTextAlign(22); // 文本居中对齐
            t->Draw();
        }
    }
    canvas_data->SaveAs("nsig_data.png");

    for(int i=0;i<13;i++){
        TCanvas *c2 = new TCanvas(Form("c2_%d",i),Form("P_bin_%d",i+1),1200,900);
        
        TGraphErrors* gr_eff = new TGraphErrors();
        int pointCount = 0;
        
        for(int j=0;j<10;j++){
            if(nsig[i][j] > 0 && nsigError[i][j] > 0){
                double x = -0.9+0.2*j; // cosθ* bin中心
                double y = nsig[i][j];
                double ex = 0.1; // bin宽度的一半
                double ey = nsigError[i][j];
                
                gr_eff->SetPoint(pointCount, x, y);
                gr_eff->SetPointError(pointCount, ex, ey);
                pointCount++;
            }
        }
        
        if(pointCount > 0){
            gr_eff->SetTitle(Form("Nsig(data) %.1f<P()Ks<%.1f;cos#theta*;Events", (i+1)/10.0, (i+2)/10.0));
            gr_eff->SetMarkerStyle(20);
            gr_eff->SetMarkerColor(kRed);
            gr_eff->SetMarkerSize(1.2);
            gr_eff->SetLineColor(kRed);
            gr_eff->Draw("AP");

            TLegend* leg = new TLegend(0.8,0.8,0.9,0.9);
            leg->SetBorderSize(1);
            leg->SetFillColor(0);
            leg->AddEntry(gr_eff, "Data", "lep");
            leg->SetTextSize(0.03);
            leg->Draw();
            
            // 设置坐标轴范围
            gr_eff->GetXaxis()->SetRangeUser(-1, 1);
            gr_eff->GetYaxis()->SetRangeUser(0, 1.2*max1[i]);
            
            c2->SaveAs(Form("nsig%d.png", i+1));
        }
        
        delete gr_eff;
        delete c1;
    }

    std::ofstream rho00file("rho00_results.txt");
    if (rho00file.is_open()) {
        rho00file << "Momentum_Bin Momentum_Center(GeV/c) rho00 rho00_error" << std::endl;
        for (int i = 0; i < 13; i++) {
            rho00file << (i + 1) << " " << momentum_center[i] << " " 
                     << rho00_values[i] << " " << rho00_errors[i] << std::endl;
        }
        rho00file.close();
        cout << "rho00_results.txt written successfully." << endl;
    } else {
        cerr << "Error: Cannot open rho00_results.txt" << endl;
    }
    
    // 4. 绘制rho00随动量变化的点图
    TCanvas* rho00_canvas = new TCanvas("rho00_canvas", "rho00 vs Momentum", 800, 600);
    TGraphErrors* rho00_graph = new TGraphErrors(13, momentum_center, rho00_values, 0, rho00_errors);
    
    // 设置图形样式
    rho00_graph->SetTitle("#rho_{00};P(Ks) [GeV/c];#rho_{00}");
    rho00_graph->SetMarkerStyle(20);
    rho00_graph->SetMarkerColor(kBlue);
    rho00_graph->SetMarkerSize(0.6);
    rho00_graph->SetLineColor(kBlue);
    rho00_graph->GetYaxis()->SetRangeUser(0.0, 0.6); // 设置Y轴范围
    rho00_graph->GetXaxis()->SetRangeUser(0, 1.4); // 设置X轴范围

    //TF1* fit_func = new TF1("fit_func", "[0]", 0, 1.4);
    //fit_func->SetParameter(0, 0.33);
    //rho00_graph->Fit("fit_func", "Q");
    
    // 绘制点图
    rho00_graph->Draw("AP");
    
    // 添加纵坐标为0.33的虚线（修改变量名避免冲突）
    TLine* ref_line = new TLine(0, 0.33, 1.4, 0.333);  // 将line改为ref_line
    ref_line->SetLineColor(kRed);
    ref_line->SetLineStyle(2); // 虚线样式
    ref_line->SetLineWidth(2);
    ref_line->Draw("same");
/*
    int n_points1 = 100;
    TGraphErrors* error_band1 = new TGraphErrors(n_points1);
    double x_min1 = rho00_graph->GetXaxis()->GetXmin();
    double x_max1 = rho00_graph->GetXaxis()->GetXmax();
    cout << "x_min1: " << x_min1 << ", x_max1: " << x_max1 << endl;
    for (int j = 0; j < n_points1; j++) {
        double x = x_min1 + j * (x_max1 - x_min1) / (n_points1 - 1);
        double y = fit_func->Eval(x);

        //cout << "Point " << j << ": x = " << x << endl;
        error_band1->SetPoint(j, x, y);
        error_band1->SetPointError(j ,0 ,0);
    }
    
    TVirtualFitter* fitter1 = TVirtualFitter::GetFitter();
    if (fitter1) {
        fitter1->GetConfidenceIntervals(error_band1, 0.683); // 1 sigma
    }
    
    error_band1->SetFillColorAlpha(kBlue - 7, 0.35);
    error_band1->SetLineColor(kBlue - 7);
    error_band1->Draw("3 same");
*/

    // 添加图例
    TLegend* rho00_legend = new TLegend(0.65, 0.75, 0.85, 0.85);
    rho00_legend->AddEntry(rho00_graph, "#rho_{00} Data", "lep");
    rho00_legend->AddEntry(ref_line, "#rho_{00} = 0.33", "l");
    rho00_legend->SetBorderSize(0);
    //rho00_legend->AddEntry(error_band1, "1 #sigma error band", "f");
    rho00_legend->Draw();
    
    // 保存图片
    rho00_canvas->SaveAs("rho00.png");
    
    // 清理内存
    delete rho00_legend;
    delete ref_line;  // 变量名相应修改
    delete rho00_graph;
    delete rho00_canvas;

    cout << "nsig_real_results.txt written successfully." << endl;
    cout << "nsig_real_2D_distribution.png written successfully." << endl;
}

