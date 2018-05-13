#include <windows.h>
#include <stdio.h>
#include "../backend/Lab2data.h"

LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);

int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
	PWSTR pCmdLine, int nCmdShow) {

	MSG  msg;
	WNDCLASSW wc = { 0 };
	wc.lpszClassName = L"Center";
	wc.hInstance = hInstance;
	wc.hbrBackground = GetSysColorBrush(COLOR_3DFACE);
	wc.lpfnWndProc = WndProc;
	wc.hCursor = LoadCursor(0, IDC_ARROW);

	RegisterClassW(&wc);
	CreateWindowW(wc.lpszClassName, L"lab2",
		WS_OVERLAPPEDWINDOW | WS_VISIBLE,
		100, 100, 500, 500, 0, 0, hInstance, 0);

	while (GetMessage(&msg, NULL, 0, 0)) {

		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}

	return (int)msg.wParam;
}

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg,
	WPARAM wParam, LPARAM lParam) {
	FILE *ifs, *cfs, *ofs;
	static HWND e1, e2, e3;
	TCHAR ifn[100], cfn[100], ofn[100];
	int t, s_num, rule, l_upper, screen;
	double beta, delta;
	switch (msg) {
	case WM_CREATE:
		e1 = CreateWindowW(L"edit", L"", WS_CHILD | WS_VISIBLE | WS_TABSTOP, 250, 205, 200, 25, hwnd, 0, 0, 0);
		e2 = CreateWindowW(L"edit", L"", WS_CHILD | WS_VISIBLE | WS_TABSTOP, 250, 235, 200, 25, hwnd, 0, 0, 0);
		e3 = CreateWindowW(L"edit", L"", WS_CHILD | WS_VISIBLE | WS_TABSTOP, 250, 265, 200, 25, hwnd, 0, 0, 0);
		CreateWindowW(L"static", L"Data File Name", WS_CHILD | WS_VISIBLE, 80, 205, 170, 25, hwnd, 0, 0, 0);
		CreateWindowW(L"static", L"Parameter File Name", WS_CHILD | WS_VISIBLE, 80, 235, 170, 25, hwnd, 0, 0, 0);
		CreateWindowW(L"static", L"Output File Name", WS_CHILD | WS_VISIBLE, 80, 265, 170, 25, hwnd, 0, 0, 0);
		CreateWindowW(L"button", L"RUN", WS_CHILD | WS_VISIBLE | BS_PUSHBUTTON, 225, 400, 50, 25, hwnd, 0, 0, 0);
		break;
	case WM_COMMAND:
		switch (HIWORD(wParam)) {
		case BN_CLICKED:
			GetWindowText(e1, ifn, 100);
			GetWindowText(e2, cfn, 100);
			GetWindowText(e3, ofn, 100);
			ifs = fopen(ifn, "r");
			cfs = fopen(cfn, "r");
			ofs = fopen(ofn, "w");
			fscanf(cfs, "%d\n%d\n%lf\n%d\n%lf\n%d\n%d\n", &t, &s_num, &delta, &rule,
				&beta, &l_upper, &screen);
			fclose(cfs);
			double data[S_MAX];
			batch_means(ifs, ofs, t, s_num, data, delta, rule, beta, l_upper, screen);
			fclose(ifs);
			fclose(ofs);
			break;
		}
		
		break;
	case WM_DESTROY:
		PostQuitMessage(0);
		break;
	}

	return DefWindowProcW(hwnd, msg, wParam, lParam);
}


