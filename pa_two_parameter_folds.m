function results = pa_two_parameter_folds()
% PA_TWO_PARAMETER_FOLDS Continue saddle-node folds in two parameter planes.
%
%   RESULTS = PA_TWO_PARAMETER_FOLDS() computes the healthy-saddle fold of
%   the full four-state PA model in (z_fixed, sigma) and (z_fixed, g_z).
%   The fold is defined by the augmented system
%
%       f(x,p) = 0,  J(x,p)*v = 0,  v'*v = 1,
%
%   and is followed with a pseudo-arclength predictor-corrector method.
%   Equilibrium-count maps classify healthy-only, bistable, and disease-only
%   regions. Representative points are validated using all scalar-reduction
%   roots and the eigenvalues of the full four-state Jacobian.
%
%   Files written beside this function:
%       pa_two_parameter_folds.png
%       pa_two_parameter_folds.pdf
%       pa_two_parameter_fold_data.mat
%       pa_two_parameter_validation.csv

    close all;
    outputFolder = fileparts(mfilename("fullpath"));
    baseParameters = pa_params();
    solverOptions = optimoptions("fsolve", ...
        Display="off", ...
        FunctionTolerance=1.0e-11, ...
        OptimalityTolerance=1.0e-11, ...
        StepTolerance=1.0e-11, ...
        FiniteDifferenceType="central", ...
        MaxIterations=500, ...
        MaxFunctionEvaluations=30000);

    initialFold = solveInitialFold(baseParameters, solverOptions);
    fprintf("\nInitial augmented fold solve\n");
    fprintf("  z_fixed = %.9f, sigma = %.6f, g_z = %.6f\n", ...
        initialFold(5), baseParameters.sodium, baseParameters.gz);
    reportFoldQuality(initialFold, "sigma", baseParameters);

    sigmaSettings = continuationSettings("sigma");
    sigmaFold = traceFold(initialFold, "sigma", 1, baseParameters, ...
        solverOptions, sigmaSettings);
    sigmaFold = orientNullVectors(sigmaFold);

    gzInitialFold = initialFold;
    gzInitialFold(6) = baseParameters.gz;
    gzSettings = continuationSettings("gz");
    gzFoldLower = traceFold(gzInitialFold, "gz", -1, baseParameters, ...
        solverOptions, gzSettings);
    gzFoldUpper = traceFold(gzInitialFold, "gz", 1, baseParameters, ...
        solverOptions, gzSettings);
    gzFold = [fliplr(gzFoldLower(:, 2:end)), gzFoldUpper];
    gzFold = removeNearDuplicatePoints(gzFold, gzSettings.scale);
    gzFold = orientNullVectors(gzFold);

    [cuspPoint, cuspDiagnostics] = detectAndRefineCusp( ...
        gzFold, baseParameters, solverOptions, gzSettings);

    sigmaGrid = classificationGrid("sigma", baseParameters);
    gzGrid = classificationGrid("gz", baseParameters);

    validationTable = validateRepresentativePoints( ...
        sigmaFold, gzFold, sigmaGrid, gzGrid, baseParameters);

    figureHandle = createFoldFigure( ...
        sigmaFold, gzFold, cuspPoint, sigmaGrid, gzGrid);

    pngPath = fullfile(outputFolder, "pa_two_parameter_folds.png");
    pdfPath = fullfile(outputFolder, "pa_two_parameter_folds.pdf");
    dataPath = fullfile(outputFolder, "pa_two_parameter_fold_data.mat");
    validationPath = fullfile(outputFolder, "pa_two_parameter_validation.csv");

    exportgraphics(figureHandle, pngPath, Resolution=300);
    exportgraphics(figureHandle, pdfPath, ContentType="vector");
    writetable(validationTable, validationPath);

    results = struct();
    results.sigmaFold = foldTable(sigmaFold, "sigma", baseParameters);
    results.gzFold = foldTable(gzFold, "gz", baseParameters);
    results.cuspPoint = cuspPoint;
    results.cuspDiagnostics = cuspDiagnostics;
    results.sigmaClassification = sigmaGrid;
    results.gzClassification = gzGrid;
    results.validation = validationTable;
    results.outputFiles = struct( ...
        png=pngPath, pdf=pdfPath, data=dataPath, validation=validationPath);
    save(dataPath, "results");

    fprintf("\nTwo-parameter continuation complete\n");
    fprintf("  sigma fold points: %d\n", size(sigmaFold, 2));
    fprintf("  g_z fold points:   %d\n", size(gzFold, 2));
    if isempty(cuspPoint)
        fprintf("  cusp: no verified cusp on the traced physical-domain fold\n");
    else
        fprintf("  cusp: z_fixed = %.8f, g_z = %.8f\n", ...
            cuspPoint(5), cuspPoint(6));
        nondegeneracy = cuspDiagnostics.nondegeneracy;
        fprintf("    scalar derivatives R_a=%.3e, R_aa=%.3e, R_aaa=%.3e\n", ...
            nondegeneracy.firstAldosteroneDerivative, ...
            nondegeneracy.secondAldosteroneDerivative, ...
            nondegeneracy.thirdAldosteroneDerivative);
        fprintf("    unfolding determinant=%.3e; next J singular value=%.3e\n", ...
            nondegeneracy.unfoldingDeterminant, ...
            nondegeneracy.jacobianSingularValues(2));
    end
    fprintf("  figure:     %s\n", pngPath);
    fprintf("  fold data:  %s\n", dataPath);
    fprintf("  validation: %s\n", validationPath);
end

function foldPoint = solveInitialFold(baseParameters, solverOptions)
    parameters = baseParameters;
    parameters.zfixed = 0.36;
    parameters.sodium = 0.0;

    healthySeed = [0.99; 0.25; 0.00; saltTarget(0.25, parameters)];
    equilibriumFunction = @(state) modelResidualAndJacobian(state, parameters);
    healthyState = fsolve(@(state) firstOutput(equilibriumFunction, state), ...
        healthySeed, solverOptions);

    [~, jacobian] = modelResidualAndJacobian(healthyState, parameters);
    [rightVectors, eigenvalues] = eig(jacobian);
    [~, nullIndex] = min(abs(diag(eigenvalues)));
    nullVector = real(rightVectors(:, nullIndex));
    nullVector = nullVector/norm(nullVector);

    fixedSigma = 0.0;
    initialUnknown = [healthyState; 0.367; nullVector];
    residualFunction = @(unknown) fixedParameterFoldResidual( ...
        unknown, fixedSigma, "sigma", baseParameters);
    [solution, residual, exitFlag] = fsolve( ...
        residualFunction, initialUnknown, solverOptions);

    if (exitFlag <= 0) || (norm(residual, inf) > 1.0e-7)
        error("pa:fold:InitialSolveFailed", ...
            "The initial fold solve did not converge. Residual infinity norm: %.3e.", ...
            norm(residual, inf));
    end

    foldPoint = [solution(1:5); fixedSigma; solution(6:9)];
end

function residual = fixedParameterFoldResidual( ...
        unknown, fixedParameter, planeName, baseParameters)
    foldPoint = [unknown(1:5); fixedParameter; unknown(6:9)];
    residual = foldResidual(foldPoint, planeName, baseParameters);
end

function settings = continuationSettings(planeName)
    settings = struct();
    settings.initialStep = 0.025;
    settings.minimumStep = 0.001;
    settings.maximumStep = 0.060;
    settings.maximumPoints = 420;
    settings.correctorTolerance = 2.0e-7;

    switch planeName
        case "sigma"
            settings.scale = [1.0; 5.0; 2.0; 1.0; 0.5; 0.4; ...
                1.0; 1.0; 1.0; 1.0];
            settings.parameterBounds = [-0.03, 0.46];
            settings.zfixedBounds = [-0.08, 0.70];
        case "gz"
            settings.scale = [1.0; 5.0; 4.0; 1.0; 0.8; 4.0; ...
                1.0; 1.0; 1.0; 1.0];
            settings.parameterBounds = [0.03, 8.0];
            settings.zfixedBounds = [-0.15, 1.50];
        otherwise
            error("pa:fold:UnknownPlane", ...
                "Unknown continuation plane '%s'.", planeName);
    end
end

function points = traceFold( ...
        initialPoint, planeName, direction, baseParameters, ...
        solverOptions, settings)
    scale = settings.scale;
    scaledResidual = @(scaledPoint) foldResidual( ...
        scale.*scaledPoint, planeName, baseParameters);

    currentPoint = initialPoint./scale;
    currentTangent = curveTangent(scaledResidual, currentPoint);
    if direction*currentTangent(6) < 0.0
        currentTangent = -currentTangent;
    end

    points = nan(10, settings.maximumPoints);
    points(:, 1) = initialPoint;
    pointCount = 1;
    stepSize = settings.initialStep;

    while pointCount < settings.maximumPoints
        predictedPoint = currentPoint + stepSize*currentTangent;
        corrector = @(candidate) [ ...
            scaledResidual(candidate); ...
            currentTangent'*(candidate - predictedPoint)];
        [correctedPoint, residual, exitFlag, output] = fsolve( ...
            corrector, predictedPoint, solverOptions);

        isAccepted = (exitFlag > 0) && ...
            (norm(residual, inf) <= settings.correctorTolerance);
        if ~isAccepted
            stepSize = 0.5*stepSize;
            if stepSize < settings.minimumStep
                fprintf("  %s continuation stopped: corrector failed at %d points.\n", ...
                    planeName, pointCount);
                break;
            end
            continue;
        end

        physicalPoint = scale.*correctedPoint;
        if physicalPoint(7:10)'*points(7:10, pointCount) < 0.0
            physicalPoint(7:10) = -physicalPoint(7:10);
            correctedPoint(7:10) = -correctedPoint(7:10);
        end

        if ~isInsideContinuationBounds(physicalPoint, settings)
            break;
        end

        newTangent = curveTangent(scaledResidual, correctedPoint);
        if newTangent'*currentTangent < 0.0
            newTangent = -newTangent;
        end

        pointCount = pointCount + 1;
        points(:, pointCount) = physicalPoint;
        currentPoint = correctedPoint;
        currentTangent = newTangent;

        if output.iterations <= 4
            stepSize = min(1.25*stepSize, settings.maximumStep);
        elseif output.iterations >= 10
            stepSize = max(0.65*stepSize, settings.minimumStep);
        end

        if hasClosedCurve(points(:, 1:pointCount), scale)
            fprintf("  %s continuation closed after %d points.\n", ...
                planeName, pointCount);
            break;
        end
    end

    points = points(:, 1:pointCount);
end

function isInside = isInsideContinuationBounds(point, settings)
    parameterValue = point(6);
    zfixed = point(5);
    states = point(1:4);

    isInside = (parameterValue >= settings.parameterBounds(1)) && ...
        (parameterValue <= settings.parameterBounds(2)) && ...
        (zfixed >= settings.zfixedBounds(1)) && ...
        (zfixed <= settings.zfixedBounds(2)) && ...
        all(isfinite(point)) && ...
        all(states > -0.25) && ...
        all(states < 25.0);
end

function isClosed = hasClosedCurve(points, scale)
    minimumLoopPoints = 60;
    if size(points, 2) < minimumLoopPoints
        isClosed = false;
        return;
    end

    scaledDistance = norm((points(:, end) - points(:, 1))./scale);
    isClosed = scaledDistance < 0.035;
end

function tangent = curveTangent(residualFunction, point)
    residualJacobian = numericalJacobian(residualFunction, point);
    [~, ~, rightVectors] = svd(residualJacobian);
    tangent = rightVectors(:, end);
    tangent = tangent/norm(tangent);
end

function residual = foldResidual(foldPoint, planeName, baseParameters)
    state = foldPoint(1:4);
    zfixed = foldPoint(5);
    secondParameter = foldPoint(6);
    nullVector = foldPoint(7:10);

    parameters = parametersAtPoint( ...
        baseParameters, planeName, zfixed, secondParameter);
    [modelResidual, jacobian] = modelResidualAndJacobian(state, parameters);
    rowScale = [1.0; 1.0; parameters.eps; 1.0/parameters.tau_s];

    residual = [ ...
        modelResidual./rowScale; ...
        (jacobian*nullVector)./rowScale; ...
        nullVector'*nullVector - 1.0];
end

function parameters = parametersAtPoint( ...
        baseParameters, planeName, zfixed, secondParameter)
    parameters = baseParameters;
    parameters.zfixed = zfixed;
    switch planeName
        case "sigma"
            parameters.sodium = secondParameter;
        case "gz"
            parameters.gz = secondParameter;
        otherwise
            error("pa:fold:UnknownPlane", ...
                "Unknown continuation plane '%s'.", planeName);
    end
end

function [residual, jacobian] = modelResidualAndJacobian(state, parameters)
    r = state(1);
    a = state(2);
    zp = state(3);
    s = state(4);

    salt = saltTarget(a, parameters);
    activeFraction = parameters.fmin + (1.0 - parameters.fmin)* ...
        hillValue(s, parameters.Kf, parameters.nf);
    sustainFraction = hillValue(s, parameters.Kz, parameters.nz);
    reninFraction = hillValue(r, parameters.Kr, parameters.pr);

    residual = [ ...
        parameters.r0*(1.0 - s) - parameters.kr*r; ...
        parameters.beta*reninFraction + ...
            parameters.alpha*(parameters.zfixed + zp*activeFraction) - ...
            parameters.ka*a; ...
        parameters.eps*(-parameters.kz*zp + parameters.gz*sustainFraction); ...
        (salt - s)/parameters.tau_s];

    activeDerivative = (1.0 - parameters.fmin)* ...
        hillDerivative(s, parameters.Kf, parameters.nf);
    sustainDerivative = hillDerivative(s, parameters.Kz, parameters.nz);
    reninDerivative = hillDerivative(r, parameters.Kr, parameters.pr);
    saltDerivative = saltTargetDerivative(a, parameters);

    jacobian = zeros(4, 4);
    jacobian(1, 1) = -parameters.kr;
    jacobian(1, 4) = -parameters.r0;
    jacobian(2, 1) = parameters.beta*reninDerivative;
    jacobian(2, 2) = -parameters.ka;
    jacobian(2, 3) = parameters.alpha*activeFraction;
    jacobian(2, 4) = parameters.alpha*zp*activeDerivative;
    jacobian(3, 3) = -parameters.eps*parameters.kz;
    jacobian(3, 4) = parameters.eps*parameters.gz*sustainDerivative;
    jacobian(4, 2) = saltDerivative/parameters.tau_s;
    jacobian(4, 4) = -1.0/parameters.tau_s;
end

function value = hillValue(input, halfSaturation, exponent)
    inputPower = input.^exponent;
    halfPower = halfSaturation^exponent;
    value = inputPower./(halfPower + inputPower);
end

function derivative = hillDerivative(input, halfSaturation, exponent)
    halfPower = halfSaturation^exponent;
    denominator = halfPower + input.^exponent;
    derivative = exponent*halfPower*input.^(exponent - 1)./denominator.^2;
end

function target = saltTarget(aldosterone, parameters)
    rawTarget = hillValue(aldosterone, parameters.Ks, parameters.ns) + ...
        parameters.sodium;
    target = min(max(rawTarget, 0.0), 0.999);
end

function derivative = saltTargetDerivative(aldosterone, parameters)
    rawTarget = hillValue(aldosterone, parameters.Ks, parameters.ns) + ...
        parameters.sodium;
    if (rawTarget <= 0.0) || (rawTarget >= 0.999)
        derivative = 0.0;
    else
        derivative = hillDerivative( ...
            aldosterone, parameters.Ks, parameters.ns);
    end
end

function jacobian = numericalJacobian(functionHandle, point)
    functionValue = functionHandle(point);
    rowCount = numel(functionValue);
    columnCount = numel(point);
    jacobian = zeros(rowCount, columnCount);

    for column = 1:columnCount
        step = 2.0e-5*max(1.0, abs(point(column)));
        forwardPoint = point;
        backwardPoint = point;
        forwardPoint(column) = forwardPoint(column) + step;
        backwardPoint(column) = backwardPoint(column) - step;
        jacobian(:, column) = ...
            (functionHandle(forwardPoint) - functionHandle(backwardPoint))/(2.0*step);
    end
end

function points = orientNullVectors(points)
    for index = 2:size(points, 2)
        if points(7:10, index)'*points(7:10, index - 1) < 0.0
            points(7:10, index) = -points(7:10, index);
        end
    end
end

function points = removeNearDuplicatePoints(points, scale)
    if size(points, 2) < 2
        return;
    end

    keep = true(1, size(points, 2));
    for index = 2:size(points, 2)
        scaledDistance = norm((points(:, index) - points(:, index - 1))./scale);
        keep(index) = scaledDistance > 1.0e-7;
    end
    points = points(:, keep);
end

function [cuspPoint, diagnostics] = detectAndRefineCusp( ...
        foldPoints, baseParameters, solverOptions, settings)
    coefficient = nan(1, size(foldPoints, 2));
    for index = 1:size(foldPoints, 2)
        coefficient(index) = cuspCoefficient( ...
            foldPoints(:, index), "gz", baseParameters);
    end

    finiteCoefficient = coefficient(isfinite(coefficient));
    diagnostics = struct();
    diagnostics.quadraticCoefficient = coefficient;
    diagnostics.signChangeIndices = find( ...
        coefficient(1:end-1).*coefficient(2:end) < 0.0);
    diagnostics.minimumAbsoluteCoefficient = ...
        min(abs(finiteCoefficient), [], "omitnan");
    cuspPoint = [];

    candidateIndices = diagnostics.signChangeIndices;
    if isempty(candidateIndices)
        [~, candidateIndex] = min(abs(coefficient));
        referenceScale = median(abs(finiteCoefficient), "omitnan");
        if abs(coefficient(candidateIndex)) < 0.02*referenceScale
            candidateIndices = candidateIndex;
        end
    end

    if isempty(candidateIndices)
        return;
    end

    candidateIndex = candidateIndices(1);
    if candidateIndex < size(foldPoints, 2)
        seedPoint = 0.5*( ...
            foldPoints(:, candidateIndex) + foldPoints(:, candidateIndex + 1));
    else
        seedPoint = foldPoints(:, candidateIndex);
    end

    coefficientScale = max(median(abs(finiteCoefficient), "omitnan"), 1.0e-4);
    scaledSeed = seedPoint./settings.scale;
    cuspResidual = @(scaledPoint) [ ...
        foldResidual(settings.scale.*scaledPoint, "gz", baseParameters); ...
        cuspCoefficient(settings.scale.*scaledPoint, "gz", baseParameters)/ ...
            coefficientScale];
    [scaledSolution, residual, exitFlag] = fsolve( ...
        cuspResidual, scaledSeed, solverOptions);
    candidatePoint = settings.scale.*scaledSolution;

    diagnostics.refinementExitFlag = exitFlag;
    diagnostics.refinementResidual = norm(residual, inf);
    diagnostics.refinedCoefficient = cuspCoefficient( ...
        candidatePoint, "gz", baseParameters);

    if (exitFlag > 0) && ...
            (norm(residual, inf) < 1.0e-6) && ...
            isInsideContinuationBounds(candidatePoint, settings)
        nondegeneracy = verifyCuspNondegeneracy( ...
            candidatePoint, baseParameters);
        diagnostics.nondegeneracy = nondegeneracy;
        if nondegeneracy.isVerified
            cuspPoint = candidatePoint;
            diagnostics.verified = true;
        else
            diagnostics.verified = false;
        end
    else
        diagnostics.verified = false;
    end
end

function coefficient = cuspCoefficient(foldPoint, planeName, baseParameters)
    state = foldPoint(1:4);
    nullVector = foldPoint(7:10);
    parameters = parametersAtPoint( ...
        baseParameters, planeName, foldPoint(5), foldPoint(6));
    [~, jacobian] = modelResidualAndJacobian(state, parameters);

    [leftVectors, ~, ~] = svd(jacobian);
    leftNull = leftVectors(:, end);
    normalization = leftNull'*nullVector;
    if abs(normalization) < 1.0e-10
        coefficient = NaN;
        return;
    end
    leftNull = leftNull/normalization;

    step = 2.0e-4/max(norm(nullVector), 1.0);
    forwardResidual = modelResidualAndJacobian( ...
        state + step*nullVector, parameters);
    centerResidual = modelResidualAndJacobian(state, parameters);
    backwardResidual = modelResidualAndJacobian( ...
        state - step*nullVector, parameters);
    secondDirectionalDerivative = ...
        (forwardResidual - 2.0*centerResidual + backwardResidual)/step^2;
    coefficient = 0.5*leftNull'*secondDirectionalDerivative;
end

function diagnostics = verifyCuspNondegeneracy(cuspPoint, baseParameters)
    parameters = parametersAtPoint( ...
        baseParameters, "gz", cuspPoint(5), cuspPoint(6));
    aldosterone = cuspPoint(2);
    aldosteroneScale = max(1.0, abs(aldosterone));
    firstStep = 1.0e-5*aldosteroneScale;
    secondStep = 1.5e-4*aldosteroneScale;
    thirdStep = 7.5e-4*aldosteroneScale;
    center = scalarEquilibriumResidual(aldosterone, parameters);
    firstPlus = scalarEquilibriumResidual( ...
        aldosterone + firstStep, parameters);
    firstMinus = scalarEquilibriumResidual( ...
        aldosterone - firstStep, parameters);
    firstDerivative = (firstPlus - firstMinus)/(2.0*firstStep);

    secondPlus = scalarEquilibriumResidual( ...
        aldosterone + secondStep, parameters);
    secondMinus = scalarEquilibriumResidual( ...
        aldosterone - secondStep, parameters);
    secondDerivative = ...
        (secondPlus - 2.0*center + secondMinus)/secondStep^2;

    thirdPlusOne = scalarEquilibriumResidual( ...
        aldosterone + thirdStep, parameters);
    thirdMinusOne = scalarEquilibriumResidual( ...
        aldosterone - thirdStep, parameters);
    thirdPlusTwo = scalarEquilibriumResidual( ...
        aldosterone + 2.0*thirdStep, parameters);
    thirdMinusTwo = scalarEquilibriumResidual( ...
        aldosterone - 2.0*thirdStep, parameters);
    thirdDerivative = ...
        (thirdPlusTwo - 2.0*thirdPlusOne + ...
        2.0*thirdMinusOne - thirdMinusTwo)/(2.0*thirdStep^3);

    zfixedStep = 1.0e-4*max(1.0, abs(parameters.zfixed));
    gzStep = 1.0e-4*max(1.0, abs(parameters.gz));
    mixedStep = 1.5e-4*aldosteroneScale;
    [zfixedDerivative, aldosteroneZfixedDerivative] = parameterDerivatives( ...
        aldosterone, mixedStep, parameters, "zfixed", zfixedStep);
    [gzDerivative, aldosteroneGzDerivative] = parameterDerivatives( ...
        aldosterone, mixedStep, parameters, "gz", gzStep);
    unfoldingMatrix = [ ...
        zfixedDerivative, gzDerivative; ...
        aldosteroneZfixedDerivative, aldosteroneGzDerivative];

    [~, jacobian] = modelResidualAndJacobian( ...
        cuspPoint(1:4), parameters);
    singularValues = sort(svd(jacobian));

    diagnostics = struct();
    diagnostics.scalarResidual = center;
    diagnostics.firstAldosteroneDerivative = firstDerivative;
    diagnostics.secondAldosteroneDerivative = secondDerivative;
    diagnostics.thirdAldosteroneDerivative = thirdDerivative;
    diagnostics.unfoldingMatrix = unfoldingMatrix;
    diagnostics.unfoldingDeterminant = det(unfoldingMatrix);
    diagnostics.jacobianSingularValues = singularValues;
    diagnostics.isVerified = ...
        (abs(center) < 1.0e-7) && ...
        (abs(firstDerivative) < 1.0e-6) && ...
        (abs(secondDerivative) < 1.0e-4) && ...
        (abs(thirdDerivative) > 1.0e-3) && ...
        (abs(diagnostics.unfoldingDeterminant) > 1.0e-4) && ...
        (singularValues(1) < 1.0e-7) && ...
        (singularValues(2) > 1.0e-3);
end

function [parameterDerivative, mixedDerivative] = parameterDerivatives( ...
        aldosterone, aldosteroneStep, parameters, fieldName, parameterStep)
    plusParameters = parameters;
    minusParameters = parameters;
    plusParameters.(fieldName) = parameters.(fieldName) + parameterStep;
    minusParameters.(fieldName) = parameters.(fieldName) - parameterStep;

    plusResidual = scalarEquilibriumResidual(aldosterone, plusParameters);
    minusResidual = scalarEquilibriumResidual(aldosterone, minusParameters);
    parameterDerivative = ...
        (plusResidual - minusResidual)/(2.0*parameterStep);

    derivativeAtPlusA = ( ...
        scalarEquilibriumResidual( ...
            aldosterone + aldosteroneStep, plusParameters) - ...
        scalarEquilibriumResidual( ...
            aldosterone + aldosteroneStep, minusParameters))/ ...
        (2.0*parameterStep);
    derivativeAtMinusA = ( ...
        scalarEquilibriumResidual( ...
            aldosterone - aldosteroneStep, plusParameters) - ...
        scalarEquilibriumResidual( ...
            aldosterone - aldosteroneStep, minusParameters))/ ...
        (2.0*parameterStep);
    mixedDerivative = ...
        (derivativeAtPlusA - derivativeAtMinusA)/(2.0*aldosteroneStep);
end

function grid = classificationGrid(planeName, baseParameters)
    switch planeName
        case "sigma"
            parameterValues = linspace(0.0, 0.40, 101);
            zfixedValues = linspace(0.0, 0.55, 111);
        case "gz"
            parameterValues = linspace(0.10, 6.0, 119);
            zfixedValues = linspace(0.0, 0.90, 109);
        otherwise
            error("pa:fold:UnknownPlane", ...
                "Unknown classification plane '%s'.", planeName);
    end

    aldosteroneGrid = linspace(0.0, 25.0, 2001);
    regionClass = zeros(numel(zfixedValues), numel(parameterValues), "uint8");
    rootCount = zeros(size(regionClass), "uint8");

    fprintf("\nClassifying %s parameter plane (%d points) ...\n", ...
        planeName, numel(regionClass));
    for parameterIndex = 1:numel(parameterValues)
        for zfixedIndex = 1:numel(zfixedValues)
            parameters = parametersAtPoint(baseParameters, planeName, ...
                zfixedValues(zfixedIndex), parameterValues(parameterIndex));
            residual = scalarEquilibriumResidual(aldosteroneGrid, parameters);
            crossingIndices = find(residual(1:end-1).*residual(2:end) < 0.0);
            count = numel(crossingIndices);
            rootCount(zfixedIndex, parameterIndex) = uint8(count);

            if count >= 3
                regionClass(zfixedIndex, parameterIndex) = uint8(2);
            elseif count == 1
                approximateRoot = mean(aldosteroneGrid( ...
                    crossingIndices(1):crossingIndices(1) + 1));
                if approximateRoot < 2.0
                    regionClass(zfixedIndex, parameterIndex) = uint8(1);
                else
                    regionClass(zfixedIndex, parameterIndex) = uint8(3);
                end
            else
                regionClass(zfixedIndex, parameterIndex) = uint8(0);
            end
        end
    end

    grid = struct();
    grid.planeName = planeName;
    grid.parameterValues = parameterValues;
    grid.zfixedValues = zfixedValues;
    grid.regionClass = regionClass;
    grid.rootCount = rootCount;
end

function residual = scalarEquilibriumResidual(aldosterone, parameters)
    salt = min(max( ...
        hillValue(aldosterone, parameters.Ks, parameters.ns) + ...
        parameters.sodium, 0.0), 0.999);
    renin = parameters.r0*(1.0 - salt)/parameters.kr;
    plasticCapacity = (parameters.gz/parameters.kz)* ...
        hillValue(salt, parameters.Kz, parameters.nz);
    activeFraction = parameters.fmin + (1.0 - parameters.fmin)* ...
        hillValue(salt, parameters.Kf, parameters.nf);
    reninProduction = parameters.beta* ...
        hillValue(renin, parameters.Kr, parameters.pr);

    residual = reninProduction + ...
        parameters.alpha*(parameters.zfixed + ...
        plasticCapacity.*activeFraction) - ...
        parameters.ka*aldosterone;
end

function validationTable = validateRepresentativePoints( ...
        sigmaFold, gzFold, sigmaGrid, gzGrid, baseParameters)
    sigmaTarget = 0.10;
    [~, sigmaIndex] = min(abs(sigmaFold(6, :) - sigmaTarget));
    sigmaZfold = sigmaFold(5, sigmaIndex);
    gzTarget = baseParameters.gz;
    [~, gzIndex] = min(abs(gzFold(6, :) - gzTarget));
    gzZfold = gzFold(5, gzIndex);

    pointSpecifications = { ...
        "sigma below fold", "sigma", max(0.01, sigmaZfold - 0.05), sigmaTarget; ...
        "sigma above fold", "sigma", sigmaZfold + 0.05, sigmaTarget; ...
        "gz below z-fold", "gz", max(0.01, gzZfold - 0.05), gzTarget; ...
        "gz above z-fold", "gz", gzZfold + 0.05, gzTarget};

    pointSpecifications = [pointSpecifications; ...
        representativeClassPoint(gzGrid, 1, "gz healthy-only"); ...
        representativeClassPoint(gzGrid, 2, "gz bistable"); ...
        representativeClassPoint(gzGrid, 3, "gz disease-only")];

    maximumEquilibriaPerPoint = 5;
    rows = cell(size(pointSpecifications, 1)*maximumEquilibriaPerPoint, 10);
    rowCount = 0;
    fprintf("\nRepresentative equilibrium and stability checks\n");
    for pointIndex = 1:size(pointSpecifications, 1)
        label = pointSpecifications{pointIndex, 1};
        planeName = pointSpecifications{pointIndex, 2};
        zfixed = pointSpecifications{pointIndex, 3};
        secondParameter = pointSpecifications{pointIndex, 4};
        parameters = parametersAtPoint( ...
            baseParameters, planeName, zfixed, secondParameter);
        roots = findEquilibriumRoots(parameters);

        fprintf("\n  %s: z_fixed=%.5f, %s=%.5f, equilibria=%d\n", ...
            label, zfixed, planeName, secondParameter, numel(roots));
        for rootIndex = 1:numel(roots)
            state = stateFromAldosterone(roots(rootIndex), parameters);
            [modelResidual, jacobian] = ...
                modelResidualAndJacobian(state, parameters);
            eigenvalues = eig(jacobian);
            maximumRealEigenvalue = max(real(eigenvalues));
            isStable = maximumRealEigenvalue < -1.0e-7;
            fprintf("    a=%8.5f  max Re(lambda)=%+.4e  %s\n", ...
                roots(rootIndex), maximumRealEigenvalue, ...
                stabilityText(isStable));
            rowCount = rowCount + 1;
            if rowCount > size(rows, 1)
                error("pa:fold:UnexpectedRootCount", ...
                    "More than five equilibria were found at a validation point.");
            end
            rows(rowCount, :) = { ...
                label, planeName, zfixed, secondParameter, rootIndex, ...
                state(1), state(2), state(3), maximumRealEigenvalue, ...
                isStable};

            if norm(modelResidual, inf) > 1.0e-7
                error("pa:fold:ValidationResidual", ...
                    "Representative equilibrium residual is too large.");
            end
        end
    end

    rows = rows(1:rowCount, :);
    validationTable = cell2table(rows, VariableNames=[ ...
        "PointLabel", "Plane", "zfixed", "SecondParameter", ...
        "EquilibriumIndex", "Renin", "Aldosterone", "PlasticCapacity", ...
        "MaximumRealEigenvalue", "IsStable"]);

    validateGridAgainstFold(sigmaGrid, sigmaFold);
    validateGridAgainstFold(gzGrid, gzFold);
end

function specification = representativeClassPoint(grid, classValue, label)
    if ~any(grid.regionClass == classValue, "all")
        specification = {label, grid.planeName, NaN, NaN};
        return;
    end

    classMask = grid.regionClass == classValue;
    distance = binaryInteriorDepth(classMask);
    [~, linearIndex] = max(distance, [], "all", "linear");
    [zIndex, parameterIndex] = ind2sub(size(classMask), linearIndex);
    specification = {label, grid.planeName, ...
        grid.zfixedValues(zIndex), grid.parameterValues(parameterIndex)};
end

function depth = binaryInteriorDepth(mask)
    depth = zeros(size(mask));
    erodedMask = mask;
    iteration = 0;

    while any(erodedMask, "all")
        iteration = iteration + 1;
        depth(erodedMask) = iteration;
        up = [false(1, size(mask, 2)); erodedMask(1:end - 1, :)];
        down = [erodedMask(2:end, :); false(1, size(mask, 2))];
        left = [false(size(mask, 1), 1), erodedMask(:, 1:end - 1)];
        right = [erodedMask(:, 2:end), false(size(mask, 1), 1)];
        erodedMask = erodedMask & up & down & left & right;
    end
end

function validateGridAgainstFold(grid, foldPoints)
    physical = (foldPoints(5, :) >= min(grid.zfixedValues)) & ...
        (foldPoints(5, :) <= max(grid.zfixedValues)) & ...
        (foldPoints(6, :) >= min(grid.parameterValues)) & ...
        (foldPoints(6, :) <= max(grid.parameterValues));
    if ~any(physical)
        error("pa:fold:MissingPhysicalFold", ...
            "No fold points lie in the %s classification grid.", grid.planeName);
    end
end

function roots = findEquilibriumRoots(parameters)
    aldosteroneGrid = linspace(0.0, 25.0, 10001);
    residual = scalarEquilibriumResidual(aldosteroneGrid, parameters);
    crossingIndices = find(residual(1:end-1).*residual(2:end) < 0.0);
    roots = nan(1, numel(crossingIndices));

    for index = 1:numel(crossingIndices)
        bracket = aldosteroneGrid(crossingIndices(index):crossingIndices(index) + 1);
        roots(index) = fzero( ...
            @(a) scalarEquilibriumResidual(a, parameters), bracket);
    end
    roots = unique(round(roots, 10));
end

function state = stateFromAldosterone(aldosterone, parameters)
    salt = saltTarget(aldosterone, parameters);
    renin = parameters.r0*(1.0 - salt)/parameters.kr;
    plasticCapacity = (parameters.gz/parameters.kz)* ...
        hillValue(salt, parameters.Kz, parameters.nz);
    state = [renin; aldosterone; plasticCapacity; salt];
end

function text = stabilityText(isStable)
    if isStable
        text = "stable";
    else
        text = "unstable";
    end
end

function figureHandle = createFoldFigure( ...
        sigmaFold, gzFold, cuspPoint, sigmaGrid, gzGrid)
    figureHandle = figure( ...
        Name="Two-parameter fold continuation", ...
        Color="w", ...
        Position=[60, 60, 1240, 540]);
    layout = tiledlayout(figureHandle, 1, 2, ...
        TileSpacing="compact", Padding="compact");

    regionColors = [ ...
        0.82, 0.94, 0.82; ...
        0.82, 0.87, 0.98; ...
        0.98, 0.84, 0.82];

    sigmaAxes = nexttile(layout);
    plotClassification(sigmaAxes, sigmaGrid, regionColors);
    hold(sigmaAxes, "on");
    physicalSigma = ...
        (sigmaFold(6, :) >= min(sigmaGrid.parameterValues)) & ...
        (sigmaFold(6, :) <= max(sigmaGrid.parameterValues)) & ...
        (sigmaFold(5, :) >= min(sigmaGrid.zfixedValues)) & ...
        (sigmaFold(5, :) <= max(sigmaGrid.zfixedValues));
    plot(sigmaAxes, sigmaFold(6, physicalSigma), sigmaFold(5, physicalSigma), ...
        "k-", LineWidth=2.4, DisplayName="healthy-saddle fold");
    xlabel(sigmaAxes, "dietary sodium bias \sigma");
    ylabel(sigmaAxes, "suppression-resistant capacity z_{fixed}");
    title(sigmaAxes, "(z_{fixed}, \sigma) fold continuation");

    gzAxes = nexttile(layout);
    plotClassification(gzAxes, gzGrid, regionColors);
    hold(gzAxes, "on");
    physicalGz = ...
        (gzFold(6, :) >= min(gzGrid.parameterValues)) & ...
        (gzFold(6, :) <= max(gzGrid.parameterValues)) & ...
        (gzFold(5, :) >= min(gzGrid.zfixedValues)) & ...
        (gzFold(5, :) <= max(gzGrid.zfixedValues));
    plot(gzAxes, gzFold(6, physicalGz), gzFold(5, physicalGz), ...
        "k-", LineWidth=2.4, DisplayName="fold locus");
    if ~isempty(cuspPoint)
        plot(gzAxes, cuspPoint(6), cuspPoint(5), "ko", ...
            MarkerFaceColor=[1.0, 0.75, 0.1], MarkerSize=8, ...
            DisplayName="verified cusp");
    end
    xlabel(gzAxes, "salt-bridge strength g_z");
    ylabel(gzAxes, "suppression-resistant capacity z_{fixed}");
    title(gzAxes, "(z_{fixed}, g_z) fold continuation");

    colormap(figureHandle, regionColors);
    title(layout, ...
        "Two-parameter bifurcation structure of the four-state PA model", ...
        FontWeight="bold");
end

function plotClassification(axesHandle, mapData, regionColors)
    imagesc(axesHandle, mapData.parameterValues, mapData.zfixedValues, ...
        double(mapData.regionClass));
    set(axesHandle, YDir="normal");
    clim(axesHandle, [0.5, 3.5]);
    grid(axesHandle, "on");
    box(axesHandle, "on");
    xlim(axesHandle, [min(mapData.parameterValues), max(mapData.parameterValues)]);
    ylim(axesHandle, [min(mapData.zfixedValues), max(mapData.zfixedValues)]);

    hold(axesHandle, "on");
    healthyPatch = patch(axesHandle, NaN, NaN, regionColors(1, :), ...
        EdgeColor="none", DisplayName="healthy-only");
    bistablePatch = patch(axesHandle, NaN, NaN, regionColors(2, :), ...
        EdgeColor="none", DisplayName="bistable");
    diseasePatch = patch(axesHandle, NaN, NaN, regionColors(3, :), ...
        EdgeColor="none", DisplayName="disease-only");
    legend(axesHandle, [healthyPatch, bistablePatch, diseasePatch], ...
        Location="best", FontSize=8);
end

function tableOutput = foldTable(points, planeName, baseParameters)
    pointCount = size(points, 2);
    maximumRealEigenvalue = nan(pointCount, 1);
    secondSmallestEigenvalue = nan(pointCount, 1);
    residualNorm = nan(pointCount, 1);
    cuspQuadraticCoefficient = nan(pointCount, 1);

    for index = 1:pointCount
        parameters = parametersAtPoint( ...
            baseParameters, planeName, points(5, index), points(6, index));
        [residual, jacobian] = ...
            modelResidualAndJacobian(points(1:4, index), parameters);
        eigenvalues = eig(jacobian);
        sortedMagnitude = sort(abs(eigenvalues));
        maximumRealEigenvalue(index) = max(real(eigenvalues));
        secondSmallestEigenvalue(index) = sortedMagnitude(2);
        residualNorm(index) = norm(residual, inf);
        cuspQuadraticCoefficient(index) = cuspCoefficient( ...
            points(:, index), planeName, baseParameters);
    end

    tableOutput = table( ...
        points(5, :)', points(6, :)', points(1, :)', points(2, :)', ...
        points(3, :)', points(4, :)', residualNorm, ...
        maximumRealEigenvalue, secondSmallestEigenvalue, ...
        cuspQuadraticCoefficient, ...
        VariableNames=[ ...
        "zfixed", "SecondParameter", "Renin", "Aldosterone", ...
        "PlasticCapacity", "SaltState", "ResidualInfinityNorm", ...
        "MaximumRealEigenvalue", "SecondSmallestEigenvalueMagnitude", ...
        "CuspQuadraticCoefficient"]);
end

function reportFoldQuality(foldPoint, planeName, baseParameters)
    parameters = parametersAtPoint( ...
        baseParameters, planeName, foldPoint(5), foldPoint(6));
    [residual, jacobian] = ...
        modelResidualAndJacobian(foldPoint(1:4), parameters);
    singularValues = sort(svd(jacobian));
    fprintf("  equilibrium residual: %.3e\n", norm(residual, inf));
    fprintf("  smallest singular value: %.3e\n", singularValues(1));
    fprintf("  next singular value:     %.3e\n", singularValues(2));
end

function value = firstOutput(functionHandle, input)
    [value, ~] = functionHandle(input);
end
