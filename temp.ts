/**
 * Converts a raw decimal string into the equivalent upload hex format
 * @param subModel - the instrument being converted
 * @param rawDecimal one line of comma separated raw instrument data
 * @param enabledSensors array of Sensors that are enabled. For 37 this is always temperature, conductivity, pressure
 * @returns a string of hex data matching what would result from uploading the same sample from the instrument
 * @throws an Error if the instrument type is not supported for this operation
 */
export function decimalToHex(
  subModel: SubModel,
  rawDecimal: string,
  enabledSensors: Sensors[],
) {
  let hex = "";
  switch (subModel) {
    case SubModel.SBE37IM:
    case SubModel.SBE37IMP:
    case SubModel.SBE37IMPODO:
    case SubModel.SBE37SI:
    case SubModel.SBE37SIP:
    case SubModel.SBE37SM:
    case SubModel.SBE37SMP:
    case SubModel.SBE37SMPODO:
      hex = createHexFromDecimal(rawDecimal, enabledSensors, false);
      break;
    case SubModel.HydroCAT:
    case SubModel.HydroCATODO:
      hex = createHexFromDecimal(rawDecimal, enabledSensors, true);
      break;
    case SubModel.SBE39Plus:
    case SubModel.SBE39PlusIM:
      hex = createHexFrom39Decimal(rawDecimal, enabledSensors);
      break;
    case SubModel.SeaFET2:
    case SubModel.SeapHox2:
      hex = createHexFromSeaFETDecimal(rawDecimal, subModel);
      break;
    case SubModel.pHyre:
    case SubModel.pHyreCAT:
      hex = createHexFrompHyreDecimal(rawDecimal, subModel);
      break;
    // not yet supported
    case SubModel.SBE19Plus:
    default:
      throw new Error(`Raw decimal to hex is not supported for ${subModel}`);
  }
  return hex;
}

/**
 * Parse the raw decimal output and convert it to the corresponding stored hex equivalent
 * @param rawDecimal - the raw decimal string output by the instrument in real-time
 * @param enabledSensors - the list of enabled sensors for the instrument
 * @param framePresent true if first value is the frame like HCAT03720842
 * @returns a string of hex data matching what would result from uploading the same sample from the instrument
 */
export function createHexFromDecimal(
  rawDecimal: string,
  enabledSensors: Sensors[],
  framePresent = false,
): string {
  try {
    let result = "";

    const rawDecimalSplit = rawDecimal.split(",");
    if (rawDecimalSplit.length === 1) {
      throw new Error(`Invalid 37 raw decimal formatting in ${rawDecimal}`);
    }

    let i = 0;
    if (framePresent) {
      i++; // skip first value, it's the frame
    }
    // add temperature hex
    if (enabledSensors.includes(Sensors.Temperature)) {
      result += parseInt(rawDecimalSplit[i])
        .toString(16)
        .padStart(HexLength.temperature, "0");
      i++;
    }

    // add cond hex if present
    if (enabledSensors.includes(Sensors.Conductivity)) {
      result += Math.round(parseFloat(rawDecimalSplit[i]) * 256)
        .toString(16)
        .padStart(HexLength.conductivity, "0");
      i++;
    }

    // save the pressure hex values to add to the result later. Ordering is mixed between decimal and hex formats
    let pressureHex = "";
    // get pressure hex and pressure temp comp if present
    if (enabledSensors.includes(Sensors.Pressure)) {
      pressureHex = parseInt(rawDecimalSplit[i])
        .toString(16)
        .padStart(HexLength.pressure, "0");
      i++;
      // add the pressure temp comp
      pressureHex += parseInt(rawDecimalSplit[i])
        .toString(16)
        .padStart(HexLength.temperatureCompensation, "0");
      i++;
    }

    // add DO hex and DO temp if present
    if (enabledSensors.includes(Sensors.SBE63)) {
      // SMP-ODO oxygen sensor is SBE63
      result += Math.round((parseFloat(rawDecimalSplit[i]) + 10) * 100000)
        .toString(16)
        .padStart(HexLength.SBE63phase, "0");
      i++;
      // add DO temperature
      result += Math.round((parseFloat(rawDecimalSplit[i]) + 1) * 1000000)
        .toString(16)
        .padStart(HexLength.SBE63temperature, "0");
      i++;
    }

    // add the pressure hex to the result now
    result += pressureHex;

    // Datetime lives at the end
    const dateTime = DateTime.fromFormat(
      rawDecimalSplit[i].trim() + rawDecimalSplit[i + 1],
      "dd MMM yyyy HH:mm:ss",
      { zone: "UTC" },
    );
    const epochSecs = dateTime.toSeconds();
    const secondsSince2000 = epochSecs - secondsBetweenEpochAnd2000;
    result += secondsSince2000.toString(16).padStart(HexLength.time, "0");

    return result.toUpperCase();
  } catch (e) {
    throw new Error(
      `Unable to parse SBE37 raw decimal data from ${rawDecimal}`,
    );
  }
}

/**
 * Parse the raw decimal output from a 39 and convert it to the corresponding stored hex equivalent
 * @param rawDecimal - the raw decimal string output by the instrument in real-time
 * @param enabledSensors - the list of enabled sensors for the instrument
 * @returns a string of hex data matching what would result from uploading the same sample from the instrument
 */
export function createHexFrom39Decimal(
  rawDecimal: string,
  enabledSensors: Sensors[],
): string {
  try {
    let result = "";

    // save the pressure hex values to add to the result later. Ordering is mixed between decimal and hex formats
    const rawDecimalSplit = rawDecimal.split(",");
    if (rawDecimalSplit.length === 1) {
      throw new Error(`Invalid 39 raw decimal formatting in ${rawDecimal}`);
    }
    const temperatureHex = parseInt(rawDecimalSplit[0])
      .toString(16)
      .padStart(HexLength.temperature, "0");

    let i = 1;
    // save the pressure hex values to add to the result later. Ordering is mixed between decimal and hex formats
    let pressureHex = "";
    // get pressure hex and pressure temp comp if present
    if (enabledSensors.includes(Sensors.Pressure)) {
      pressureHex = parseInt(rawDecimalSplit[i])
        .toString(16)
        .padStart(HexLength.pressure, "0");
      i++;
      // add the pressure temp comp
      pressureHex += parseInt(rawDecimalSplit[i])
        .toString(16)
        .padStart(HexLength.temperatureCompensation, "0");
      i++;
    }

    // convert the time to hex and start the result with it
    const dateTime = DateTime.fromFormat(
      rawDecimalSplit[i].trim() + rawDecimalSplit[i + 1],
      "dd MMM yyyy HH:mm:ss",
      { zone: "UTC" },
    );
    const epochSecs = dateTime.toSeconds();
    const secondsSince2000 = epochSecs - secondsBetweenEpochAnd2000;
    result += secondsSince2000.toString(16).padStart(HexLength.time, "0");
    result += temperatureHex;
    result += pressureHex;

    return result.toUpperCase();
  } catch (e) {
    throw new Error(
      `Unable to parse SBE39 raw decimal data from ${rawDecimal}`,
    );
  }
}

/**
 * Parse the raw decimal output from a SeaFET and convert it to the corresponding stored hex equivalent
 * @param rawDecimal - the raw decimal string output by the instrument in real-time
 * @param {SubModel} subModel - the subModel of instrument that this is (SeaFET or SeapHOx)
 * @param {boolean} [isShallow=true] - true if it's a shallow SeapHOx or SeaFET, false if it's a deep SeapHOX
 * @returns a string of hex data matching what would result from uploading the same sample from the instrument
 */
export function createHexFromSeaFETDecimal(
  rawDecimal: string,
  subModel: SubModel,
  isShallow = true,
): string {
  let result = "";

  if (
    subModel !== SubModel.SeaFET2 &&
    subModel &&
    subModel !== SubModel.SeapHox2
  ) {
    throw new Error(`${subModel.toString()} is not recognized as a SeaFET`);
  }

  try {
    const parsedRawDecimal = readRawDecimal(subModel, rawDecimal, []);

    if (subModel === SubModel.SeapHox2) {
      // add temperature
      result += parsedRawDecimal[HexDataTypes.temperature]
        .toString(16)
        .padStart(HexLength.temperature, "0");

      // add conductivity
      const rawCond = Math.round(
        (parsedRawDecimal[HexDataTypes.conductivity] as number) * 256,
      );
      result += rawCond.toString(16).padStart(HexLength.conductivity, "0");

      // add pressure
      result += parsedRawDecimal[HexDataTypes.pressure]
        .toString(16)
        .padStart(HexLength.pressure, "0");

      // add temperatureCompensation
      result += parsedRawDecimal[HexDataTypes.temperatureCompensation]
        .toString(16)
        .padStart(HexLength.temperatureCompensation, "0");

      // add SBE63phase
      const raw63Phase = Math.round(
        ((parsedRawDecimal[HexDataTypes.SBE63phase] as number) + 10) * 100000,
      );
      result += raw63Phase.toString(16).padStart(HexLength.SBE63phase, "0");

      // add SBE63temperature
      const raw63Temp = Math.round(
        ((parsedRawDecimal[HexDataTypes.SBE63temperature] as number) + 1) *
          1000000,
      );
      result += raw63Temp
        .toString(16)
        .padStart(HexLength.SBE63temperature, "0");
    }

    // external pH
    result += parsedRawDecimal[HexDataTypes.vrsExternal]
      .toString(16)
      .padStart(HexLength.vrsExternal, "0");

    if (isShallow) {
      // internal pH
      result += parsedRawDecimal[HexDataTypes.vrsInternal]
        .toString(16)
        .padStart(HexLength.vrsInternal, "0");

      // pH reference temperature
      result += parsedRawDecimal[HexDataTypes.pHtemperature]
        .toString(16)
        .padStart(HexLength.pHtemperature, "0");
    }

    result += parsedRawDecimal[HexDataTypes.vk]
      .toString(16)
      .padStart(HexLength.vk, "0");

    result += parsedRawDecimal[HexDataTypes.ib]
      .toString(16)
      .padStart(HexLength.ib, "0");

    result += parsedRawDecimal[HexDataTypes.ik]
      .toString(16)
      .padStart(HexLength.ik, "0");

    // scale internal humidity from 16 bits to 12
    result += (
      (parsedRawDecimal[HexDataTypes.relativeHumidity] as number) >>> 4
    )
      .toString(16)
      .padStart(HexLength.relativeHumidity, "0");

    // scale internal temp from 16 bits to 12
    result += (
      (parsedRawDecimal[HexDataTypes.internalTemperature] as number) >>> 4
    )
      .toString(16)
      .padStart(HexLength.internalTemperature, "0");

    // datetime
    const epochSecs = (
      parsedRawDecimal[HexDataTypes.dateTime] as DateTime
    ).toSeconds();
    const secondsSince2000 = epochSecs - secondsBetweenEpochAnd2000;
    result += secondsSince2000.toString(16).padStart(HexLength.time, "0");

    result += parsedRawDecimal[HexDataTypes.errorFlag]
      .toString(16)
      .padStart(HexLength.errorFlag, "0");
  } catch (e) {
    throw new Error(
      `Unable to parse SeaFET raw decimal data from ${rawDecimal}`,
    );
  }

  return result.toUpperCase();
}

/**
 * Parse the raw decimal output from a pHyre and convert it to the corresponding stored hex equivalent
 * @param rawDecimal - the raw decimal string output by the instrument in real-time
 * @param {SubModel} subModel - the subModel of instrument that this is (pHyre or pHyreCAT)
 * @returns a string of hex data matching what would result from uploading the same sample from the instrument
 */
export function createHexFrompHyreDecimal(
  rawDecimal: string,
  subModel: SubModel,
): string {
  let result = "";

  if (subModel !== SubModel.pHyre && subModel !== SubModel.pHyreCAT) {
    throw new Error(
      `${subModel.toString()} is not recognized as a pHyre or pHyreCAT`,
    );
  }

  try {
    const parsedRawDecimal = readRawDecimal(subModel, rawDecimal, []);

    if (subModel === SubModel.pHyreCAT) {
      result += parsedRawDecimal[HexDataTypes.temperature]
        .toString(16)
        .padStart(HexLength.temperature, "0");
      result += parsedRawDecimal[HexDataTypes.conductivity]
        .toString(16)
        .padStart(HexLength.conductivity, "0");
      result += parsedRawDecimal[HexDataTypes.pressure]
        .toString(16)
        .padStart(HexLength.pressure, "0");
      result += parsedRawDecimal[HexDataTypes.temperatureCompensation]
        .toString(16)
        .padStart(HexLength.temperatureCompensation, "0");
      result += parsedRawDecimal[HexDataTypes.SBE63phase]
        .toString(16)
        .padStart(HexLength.SBE63phase, "0");
      result += parsedRawDecimal[HexDataTypes.SBE63temperature]
        .toString(16)
        .padStart(HexLength.SBE63temperature, "0");
    }

    result += parsedRawDecimal[HexDataTypes.vrsExternal]
      .toString(16)
      .padStart(HexLength.vrsExternal, "0");
    result += parsedRawDecimal[HexDataTypes.vrsInternal]
      .toString(16)
      .padStart(HexLength.vrsInternal, "0");
    result += parsedRawDecimal[HexDataTypes.pHtemperatureExternal]
      .toString(16)
      .padStart(HexLength.pHtemperature, "0");
    result += parsedRawDecimal[HexDataTypes.pHtemperatureInternal]
      .toString(16)
      .padStart(HexLength.pHtemperature, "0");
    result += parsedRawDecimal[HexDataTypes.vk]
      .toString(16)
      .padStart(HexLength.vk, "0");
    result += parsedRawDecimal[HexDataTypes.ib]
      .toString(16)
      .padStart(HexLength.ib, "0");
    result += parsedRawDecimal[HexDataTypes.ik]
      .toString(16)
      .padStart(HexLength.ik, "0");

    const epochSecs = (
      parsedRawDecimal[HexDataTypes.dateTime] as DateTime
    ).toSeconds();
    const secondsSince2000 = epochSecs - secondsBetweenEpochAnd2000;
    result += secondsSince2000.toString(16).padStart(HexLength.time, "0");
    result += parsedRawDecimal[HexDataTypes.errorFlag]
      .toString(16)
      .padStart(HexLength.errorFlag, "0");
  } catch (e) {
    throw new Error(
      `Unable to parse pHyre raw decimal data from ${rawDecimal}`,
    );
  }

  return result.toUpperCase();
}
